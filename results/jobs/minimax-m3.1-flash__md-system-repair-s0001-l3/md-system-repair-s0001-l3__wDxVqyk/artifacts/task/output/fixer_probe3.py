from pdbfixer import PDBFixer
from openmm import unit
f=PDBFixer('/task/output/protein_clean3.pdb')
print('initial',len(list(f.topology.residues())),len(list(f.topology.atoms())))
f.findMissingResidues(); print('missingRes',f.missingResidues)
f.findMissingAtoms()
for r,v in f.missingAtoms.items(): print('missingAtoms',r.id,r.name,r.insertionCode,[a.name for a in v])
f.findNonstandardResidues(); print('nonstandard',[(r.name,r.id,r.insertionCode) for r in f.nonstandardResidues])
f.addMissingAtoms(seed=2025)
print('afterHeavy',len(list(f.topology.residues())),len(list(f.topology.atoms())))
f.addMissingHydrogens(7.0)
print('afterH',len(list(f.topology.residues())),len(list(f.topology.atoms())))
for r in f.topology.residues():
 if r.name=='CYS': print('CYS',r.id,repr(r.insertionCode),[a.name for a in r.atoms()])
with open('/task/output/protein_fixed3.pdb','w') as out:
 for i,(atom,pos) in enumerate(zip(f.topology.atoms(),f.positions),1):
  r=atom.residue; el=atom.element.symbol if atom.element else '?'
  x,y,z=[float(v) for v in pos.value_in_unit(unit.angstrom)]
  out.write(f'ATOM  {i:5d} {atom.name:>4s} {r.name:>3s} {r.chain.id}{r.id:>4s}{r.insertionCode}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n')
 out.write('TER\nEND\n')
