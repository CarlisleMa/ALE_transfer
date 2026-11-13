import pickle, numpy as np
from openmm import unit, Vec3
from openmm.app import Modeller, ForceField
with open('/task/output/protein_final_state.pkl','rb') as f: top,pos=pickle.load(f)
edges=np.array([50.165,45.429,50.023])
box=unit.Quantity((Vec3(edges[0]/10,0,0),Vec3(0,edges[1]/10,0),Vec3(0,0,edges[2]/10)),unit.nanometer)
top.setPeriodicBoxVectors(box)
mod=Modeller(top,pos); ff=ForceField('amber14-all.xml')
mod.addSolvent(ff,model='tip3p',boxVectors=box,neutralize=False)
waters=[r for r in mod.topology.residues() if r.name=='HOH']
print('WATERS',len(waters),'ATOMS',len(list(mod.topology.atoms())))
print('BOX',mod.topology.getPeriodicBoxVectors())
import pickle
pickle.dump((mod.topology,mod.positions),open('/task/output/solvated_noions_state.pkl','wb'))
