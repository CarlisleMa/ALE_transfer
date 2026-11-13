from openmm.app import PDBFile, Topology, Modeller, Element
from openmm import unit, Vec3
import numpy as np, json
# Load fixed fragments.
pA=PDBFile('/task/output/fixed_A3.pdb'); pB=PDBFile('/task/output/fixed_B3.pdb')
top=Topology(); chains={}; atom_map=[]; positions=[]
for p in [pA,pB]:
    for chain in p.topology.chains():
        if chain.id not in chains: chains[chain.id]=top.addChain(chain.id)
        for res in chain.residues():
            nr=top.addResidue(res.name,chains[chain.id],res.id,res.insertionCode)
            local=[]
            for a in res.atoms():
                na=top.addAtom(a.name,a.element,nr)
                local.append((na,a.index)); positions.append(p.positions[a.index])
            atom_map.append((res,local,nr))
# Copy bonds.
bond_lookup={}
for residue,local,nr in atom_map:
    for a in residue.atoms(): bond_lookup[(residue,a.name)]=next(x[0] for x in local if p.topology.atoms()[x[1]].name==a.name)
for p in [pA,pB]:
    for a,b in p.topology.bonds():
        if (a.residue,a.name) in bond_lookup and (b.residue,b.name) in bond_lookup:
            top.addBond(bond_lookup[(a.residue,a.name)],bond_lookup[(b.residue,b.name)])
# Explicit crystallographic disulfides, using work residue IDs.
pairs=[(5,51),(14,34),(26,48)]
def residues_by_id(t): return {r.id:r for r in t.residues()}
ra=residues_by_id(pA.topology); rb=residues_by_id(pB.topology)
# map old ID to fixed topology residue object
fixed_by_work={}
for residue,local,nr in atom_map: fixed_by_work[residue.id]=nr
# Work IDs are global sequential in the fragments; locate IDs directly.
for i,j in pairs:
    rr_i=fixed_by_work[i]; rr_j=fixed_by_work[j]
    def atom(nr,name): return next(a for a in nr.atoms() if a.name==name)
    top.addBond(atom(rr_i,'SG'),atom(rr_j,'SG'))
    print('DISULFIDE',i,j)
mod=Modeller(top,positions)
# Remove any pre-existing thiol H; Modeller chooses CYX from SG-SG bonds and omits HG.
remove=[a for a in mod.topology.atoms() if a.name=='HG' and a.residue.name=='CYS']
if remove: mod.delete(remove)
# Add any hydrogens not already present and select CYX automatically.
ff=app_forcefield=None
from openmm.app import ForceField
ff=ForceField('amber14-all.xml')
# Existing PDBFixer hydrogens are retained; adding hydrogens with variants lets Modeller adjust.
mod.addHydrogens(ff,variants=[('CYX',[26])] if False else None)
# Ensure CYX matching for all six bonded cysteines by explicit variants indices.
indices=[r.index for r in mod.topology.residues() if r.name=='CYS']
variants=[('CYX',[i]) for i in indices]
# Modeller accepts one tuple per residue; set all non-CYS None.
allvars=[]
for r in mod.topology.residues(): allvars.append(('CYX',[r.index]) if r.name=='CYS' else None)
# addHydrogens already called, so directly rename bonded CYS residues to CYX for matching.
for r in mod.topology.residues():
    if r.name=='CYS':
        sg=next(a for a in r.atoms() if a.name=='SG')
        bonded=any((a.name=='SG' and b.residue is not r) or (b.name=='SG' and a.residue is not r) for a,b in mod.topology.bonds())
        if bonded: r.name='CYX'
with open('/task/output/protein_combined.pdb','w') as f:
    from openmm.app import PDBFile as PF
    PF.writeFile(mod.topology,mod.positions,f,keepIds=True)
print('COMBINED',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
print('SGBONDS',[(a.residue.id,a.residue.name,b.residue.id,b.residue.name) for a,b in mod.topology.bonds() if a.name=='SG' or b.name=='SG'])
# Save mapping work ID -> original label.
mp=json.load(open('/task/output/residue_mapping3.json'))
json.dump(mp,open('/task/output/residue_mapping_final.json','w'),indent=2)
