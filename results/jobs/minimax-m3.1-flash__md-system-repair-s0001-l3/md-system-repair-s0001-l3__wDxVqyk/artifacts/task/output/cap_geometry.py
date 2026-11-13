# Analyze the internal-break geometry to determine cap placement vectors.
from openmm.app import PDBFile
import numpy as np
p=PDBFile('/task/output/fixed_A3.pdb'); ra=list(p.topology.residues())[-1]; pos=p.positions
for name in ['N','CA','C','O']:
 a=next(a for a in ra.atoms() if a.name==name); print('A',name,a.index,np.array(pos[a.index].value_in_unit(__import__('openmm').unit.angstrom)))
p=PDBFile('/task/output/fixed_B3.pdb'); rb=list(p.topology.residues())[0]; pos=p.positions
for name in ['N','CA','C','O']:
 a=next(a for a in rb.atoms() if a.name==name); print('B',name,a.index,np.array(pos[a.index].value_in_unit(__import__('openmm').unit.angstrom)))
