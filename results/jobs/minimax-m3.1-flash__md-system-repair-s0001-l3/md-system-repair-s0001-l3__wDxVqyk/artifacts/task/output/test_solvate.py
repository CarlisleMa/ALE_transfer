import pickle, numpy as np
from openmm import unit, Vec3
from openmm.app import Modeller, ForceField
with open('/task/output/protein_final_state.pkl','rb') as f: top,pos=pickle.load(f)
edges=np.array([50.165,45.429,50.023])
box=unit.Quantity((Vec3(edges[0]/10,0,0),Vec3(0,edges[1]/10,0),Vec3(0,0,edges[2]/10)),unit.nanometer)
top.setPeriodicBoxVectors(box)
mod=Modeller(top,pos); ff=ForceField('amber14-all.xml')
volume=edges.prod(); n=int(round(volume*0.0334))
print('VOLUME',volume,'REQUESTED_WATERS',n)
mod.addSolvent(ff,model='tip3p',numAdded=n,neutralize=False)
print('AFTER',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
waters=[r for r in mod.topology.residues() if r.name=='HOH']
print('WATERS',len(waters),'IONS',[(r.name,r.id) for r in mod.topology.residues() if r.name not in ('HOH','ACE','NME') and r.name not in ['ALA']])
print('BOX',mod.topology.getPeriodicBoxVectors())
