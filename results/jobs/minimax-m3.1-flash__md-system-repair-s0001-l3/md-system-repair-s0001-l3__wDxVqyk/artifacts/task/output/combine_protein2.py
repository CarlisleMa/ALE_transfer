from openmm.app import PDBFile, Topology, Modeller, PDBFile as Writer
from openmm import unit
import json
pA=PDBFile('/task/output/fixed_A3.pdb'); pB=PDBFile('/task/output/fixed_B3.pdb')
top=Topology(); chain_map={}; atom_map={}; positions=[]; residue_map={}
for p in [pA,pB]:
    for chain in p.topology.chains():
        if chain.id not in chain_map: chain_map[chain.id]=top.addChain(chain.id)
        for oldres in chain.residues():
            newres=top.addResidue(oldres.name,chain_map[chain.id],oldres.id,oldres.insertionCode)
            residue_map[oldres]=newres
            for oldatom in oldres.atoms():
                newatom=top.addAtom(oldatom.name,oldatom.element,newres)
                atom_map[oldatom]=newatom
                positions.append(p.positions[oldatom.index])
for p in [pA,pB]:
    for a,b in p.topology.bonds(): top.addBond(atom_map[a],atom_map[b])
# Map globally unique work IDs to original PDB residue labels.
mapping=json.load(open('/task/output/residue_mapping3.json'))
work_to_orig={x['work_num']:(x['old_num'],x['old_insertion']) for x in mapping}
# Fragment-local IDs are sequential; reconstruct each residue's global work ID from chain and local ID.
work_to_old={}
for p,offset in [(pA,0),(pB,19)]:
    for oldres in p.topology.residues(): work_to_old[int(oldres.id)+offset]=oldres
# Original pairs: 22-70, 31-55, 47-67A.
for wi,wj in [(5,51),(14,34),(26,47)]:
    oi,oj=work_to_orig[wi],work_to_orig[wj]
    ri=residue_map[work_to_old[wi]]; rj=residue_map[work_to_old[wj]]
    si=next(a for a in ri.atoms() if a.name=='SG'); sj=next(a for a in rj.atoms() if a.name=='SG')
    top.addBond(si,sj)
    print('DISULFIDE',oi,oj)
mod=Modeller(top,positions)
# All six cysteines are disulfide bonded: remove HG and use the CYX template.
remove=[a for a in mod.topology.atoms() if a.name=='HG' and a.residue.name=='CYS']
print('REMOVE_HG',len(remove))
if remove: mod.delete(remove)
for r in mod.topology.residues():
    if r.name=='CYS': r.name='CYX'
with open('/task/output/protein_combined.pdb','w') as f:
    Writer.writeFile(mod.topology,mod.positions,f,keepIds=True)
print('COMBINED',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
print('SG_BONDS',[(a.residue.id,a.residue.name,b.residue.id,b.residue.name) for a,b in mod.topology.bonds() if a.name=='SG' or b.name=='SG'])
