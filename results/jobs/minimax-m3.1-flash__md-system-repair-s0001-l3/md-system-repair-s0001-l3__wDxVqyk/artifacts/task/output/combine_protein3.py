from openmm.app import PDBFile, Topology, Modeller, PDBFile as Writer
import json
pA=PDBFile('/task/output/fixed_A3.pdb'); pB=PDBFile('/task/output/fixed_B3.pdb')
top=Topology(); chain_map={}; atom_map={}; residue_map={}; positions=[]
for p in [pA,pB]:
 for chain in p.topology.chains():
  if chain.id not in chain_map: chain_map[chain.id]=top.addChain(chain.id)
  for oldres in chain.residues():
   newres=top.addResidue(oldres.name,chain_map[chain.id],oldres.id,oldres.insertionCode)
   residue_map[(chain.id,oldres.id)]=newres
   for oldatom in oldres.atoms():
    newatom=top.addAtom(oldatom.name,oldatom.element,newres)
    atom_map[oldatom]=newatom; positions.append(p.positions[oldatom.index])
for p in [pA,pB]:
 for a,b in p.topology.bonds(): top.addBond(atom_map[a],atom_map[b])
# Global work IDs: fragment A 1-19, fragment B 20-54.
pairs=[(('A','5'),('B','51')),(('A','14'),('B','34')),(('B','26'),('B','47'))]
for keyi,keyj in pairs:
 ri=residue_map[keyi]; rj=residue_map[keyj]
 si=next(a for a in ri.atoms() if a.name=='SG'); sj=next(a for a in rj.atoms() if a.name=='SG')
 top.addBond(si,sj); print('DISULFIDE',keyi,keyj)
mod=Modeller(top,positions)
remove=[a for a in mod.topology.atoms() if a.name=='HG' and a.residue.name=='CYS']
print('REMOVE_HG',len(remove))
if remove: mod.delete(remove)
for r in mod.topology.residues():
 if r.name=='CYS': r.name='CYX'
with open('/task/output/protein_combined.pdb','w') as f: Writer.writeFile(mod.topology,mod.positions,f,keepIds=True)
print('COMBINED',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
print('SG_BONDS',[(a.residue.id,a.residue.name,b.residue.id,b.residue.name) for a,b in mod.topology.bonds() if a.name=='SG' or b.name=='SG'])
