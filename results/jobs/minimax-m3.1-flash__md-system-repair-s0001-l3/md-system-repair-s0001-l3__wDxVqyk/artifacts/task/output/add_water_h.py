import pickle
from openmm.app import Modeller, ForceField
from openmm import unit
with open('/task/output/solvated_noions_state.pkl','rb') as f: top,pos=pickle.load(f)
mod=Modeller(top,pos); ff=ForceField('amber14-all.xml')
mod.addHydrogens(ff)
print('AFTER',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
w=[r for r in mod.topology.residues() if r.name=='HOH']
print('WATERS',len(w),'FIRST',[[a.name for a in r.atoms()] for r in w[:2]])
s=ff.createSystem(mod.topology,nonbondedMethod=__import__('openmm').app.PME,nonbondedCutoff=8.5*unit.angstrom)
print('PARAM_OK',s.getNumParticles())
pickle.dump((mod.topology,mod.positions),open('/task/output/solvated_h_state.pkl','wb'))
