from openmm.app import PDBFile, Topology, Modeller, PDBFile as Writer
from openmm import unit
import numpy as np
pA=PDBFile('/task/output/fixed_A3.pdb'); pB=PDBFile('/task/output/fixed_B3.pdb')
top=Topology(); chain_map={}; atom_map={}; residue_map={}; pos=[]
for p in [pA,pB]:
 for chain in p.topology.chains():
  if chain.id not in chain_map: chain_map[chain.id]=top.addChain(chain.id)
  for oldres in chain.residues():
   newres=top.addResidue(oldres.name,chain_map[chain.id],oldres.id,oldres.insertionCode)
   residue_map[(chain.id,oldres.id)]=newres
   for oldatom in oldres.atoms():
    # All six cysteines form disulfides; omit thiol HG.
    if oldres.name=='CYS' and oldatom.name=='HG': continue
    newatom=top.addAtom(oldatom.name,oldatom.element,newres)
    atom_map[oldatom]=newatom; pos.append(p.positions[oldatom.index])
for p in [pA,pB]:
 for a,b in p.topology.bonds():
  if a in atom_map and b in atom_map: top.addBond(atom_map[a],atom_map[b])
pairs=[(('A','5'),('B','51')),(('A','14'),('B','34')),(('B','26'),('B','47'))]
for keyi,keyj in pairs:
 ri=residue_map[keyi]; rj=residue_map[keyj]
 si=next(a for a in ri.atoms() if a.name=='SG'); sj=next(a for a in rj.atoms() if a.name=='SG')
 top.addBond(si,sj); print('DISULFIDE',keyi,keyj)
for r in top.residues():
 if r.name=='CYS': r.name='CYX'
positions=unit.Quantity(np.asarray([[float(v) for v in q.value_in_unit(unit.angstrom)] for q in pos],dtype=float)/10.0,unit.nanometer)
with open('/task/output/protein_combined.pdb','w') as f: Writer.writeFile(top,positions,f,keepIds=True)
print('COMBINED',len(list(top.residues())),len(list(top.atoms())))
print('SG_BONDS',[(a.residue.id,a.residue.name,b.residue.id,b.residue.name) for a,b in top.bonds() if a.name=='SG' or b.name=='SG'])
