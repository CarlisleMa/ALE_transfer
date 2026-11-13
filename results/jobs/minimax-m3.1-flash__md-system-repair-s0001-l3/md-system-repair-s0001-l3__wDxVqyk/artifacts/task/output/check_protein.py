import pickle, numpy as np
from openmm import unit, Vec3
from openmm.app import ForceField, PME
with open('/task/output/protein_final_state.pkl','rb') as f: top,pos=pickle.load(f)
xyz=np.array(pos.value_in_unit(unit.angstrom))
lo=xyz.min(axis=0); hi=xyz.max(axis=0); extent=hi-lo
edges=extent+16.0
print('PROTEIN',len(list(top.residues())),len(list(top.atoms())))
print('EXTENT_LO',lo,'EXTENT_HI',hi,'EXTENT',extent,'EDGES',edges)
print('CHARGES_TEMPLATE')
ff=ForceField('amber14-all.xml')
box=unit.Quantity((Vec3(edges[0]/10.0,0,0),Vec3(0,edges[1]/10.0,0),Vec3(0,0,edges[2]/10.0)),unit.nanometer)
top.setPeriodicBoxVectors(box)
try:
 s=ff.createSystem(top,nonbondedMethod=PME,nonbondedCutoff=8.5*unit.angstrom)
 print('PARAM_OK',s.getNumParticles())
except Exception as e:
 print('PARAM_FAIL',type(e).__name__,e)
