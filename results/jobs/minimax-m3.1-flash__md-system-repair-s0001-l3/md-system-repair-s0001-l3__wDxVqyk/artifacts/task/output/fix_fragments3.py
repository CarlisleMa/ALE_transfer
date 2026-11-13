from pdbfixer import PDBFixer
from openmm import unit
for label in ['A','B']:
 f=PDBFixer(f'/task/output/fragment_{label}3.pdb')
 print(label,'initial',len(list(f.topology.residues())),len(list(f.topology.atoms())))
 f.findMissingResidues(); f.findMissingAtoms(); f.findNonstandardResidues()
 for r,v in f.missingAtoms.items(): print('missing',r.id,r.name,[a.name for a in v])
 f.addMissingAtoms(seed=2025); f.addMissingHydrogens(7.0)
 print(label,'fixed',len(list(f.topology.residues())),len(list(f.topology.atoms())))
 with open(f'/task/output/fixed_{label}3.pdb','w') as out:
  for i,(a,p) in enumerate(zip(f.topology.atoms(),f.positions),1):
   r=a.residue; el=a.element.symbol if a.element else '?'
   x,y,z=[float(v) for v in p.value_in_unit(unit.angstrom)]
   out.write(f'ATOM  {i:5d} {a.name:>4s} {r.name:>3s} {r.chain.id}{r.id:>4s}{r.insertionCode}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n')
  out.write('TER\nEND\n')
