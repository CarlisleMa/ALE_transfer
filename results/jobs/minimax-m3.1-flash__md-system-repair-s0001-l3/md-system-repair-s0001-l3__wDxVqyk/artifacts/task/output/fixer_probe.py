from pdbfixer import PDBFixer
from openmm import unit
unit_angstrom=unit.angstrom
f=PDBFixer('/task/output/protein_clean.pdb')
print('initial',len(list(f.topology.residues())),len(list(f.topology.atoms())))
f.findMissingResidues(); print('missingRes',f.missingResidues)
f.findMissingAtoms(); print('missingAtoms',f.missingAtoms)
f.findNonstandardResidues(); print('nonstandard',[(r.name,r.id,r.insertionCode) for r in f.nonstandardResidues])
f.addMissingAtoms(seed=2025)
print('after heavy',len(list(f.topology.residues())),len(list(f.topology.atoms())))
for a,b in f.topology.bonds():
 if a.residue is not b.residue and a.residue.chain==b.residue.chain:
  print('INTERRES',a.residue.name,a.residue.id,repr(a.residue.insertionCode),a.name,'--',b.residue.name,b.residue.id,repr(b.residue.insertionCode),b.name)
f.addMissingHydrogens(7.0)
print('after H',len(list(f.topology.residues())),len(list(f.topology.atoms())))
for r in f.topology.residues():
 if r.name in ('CYS','CYX'): print('CYS',r.id,repr(r.insertionCode),[a.name for a in r.atoms()])
print('S bonds')
for a,b in f.topology.bonds():
 if a.name=='SG' or b.name=='SG': print(a.residue.id,repr(a.residue.insertionCode),a.name,'--',b.residue.id,repr(b.residue.insertionCode),b.name)
with open('/task/output/protein_fixed.pdb','w') as out:
 for i,(atom,pos) in enumerate(zip(f.topology.atoms(),f.positions),1):
  r=atom.residue; el=atom.element.symbol if atom.element else '?'
  out.write(f'ATOM  {i:5d} {atom.name:>4s} {r.name:>3s} {r.chain.id}{r.id:>4s}{r.insertionCode}   {pos.value_in_unit(unit_angstrom):8.3f}{pos.value_in_unit(unit_angstrom):8.3f}{pos.value_in_unit(unit_angstrom):8.3f}  1.00  0.00          {el:>2s}\n')
