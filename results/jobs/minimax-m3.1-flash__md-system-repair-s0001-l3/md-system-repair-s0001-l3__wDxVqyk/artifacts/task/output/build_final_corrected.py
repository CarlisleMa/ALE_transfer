import pickle, numpy as np
from openmm import unit, Vec3
from openmm.app import Modeller, ForceField, PME
from parmed.openmm.topsystem import load_topology
FF_FILES=['amber14/protein.ff14SB.xml','amber14/tip3p.xml']
EDGES=np.array([50.165,45.429,50.023])
with open('/task/output/protein_final_state.pkl','rb') as f: top,pos=pickle.load(f)
# Rigidly center solute in the exact rectangular cell; required edge=extent+16 gives 8 A per side.
xyz=np.array([[float(x) for x in q.value_in_unit(unit.angstrom)] for q in pos])
shift=(EDGES/2.0-(xyz.min(0)+xyz.max(0))/2.0)
pos=unit.Quantity((xyz+shift)/10.0,unit.nanometer)
box=unit.Quantity((Vec3(EDGES[0]/10,0,0),Vec3(0,EDGES[1]/10,0),Vec3(0,0,EDGES[2]/10)),unit.nanometer)
top.setPeriodicBoxVectors(box)
mod=Modeller(top,pos); ff=ForceField(*FF_FILES)
mod.addSolvent(ff,model='tip3p',boxVectors=box,neutralize=False)
mod.addHydrogens(ff)
waters=[r for r in mod.topology.residues() if r.name=='HOH']
print('INITIAL_WATERS',len(waters),'BOX',mod.topology.getPeriodicBoxVectors())
# Protein is +5: 5 Cl- neutralize; 9 NaCl pairs -> 9 Na+, 14 Cl-.
water_pos={r:mod.positions[next(a for a in r.atoms() if a.name=='O').index] for r in waters}
ionic=9*55.4/(len(waters)-5)*unit.molar
mod._addIons(ff,len(waters),water_pos,positiveIon='Na+',negativeIon='Cl-',ionicStrength=ionic,neutralize=True,residueTemplates={})
final_waters=[r for r in mod.topology.residues() if r.name=='HOH']; na=[r for r in mod.topology.residues() if r.name=='NA']; cl=[r for r in mod.topology.residues() if r.name=='CL']
print('FINAL_COUNTS',len(final_waters),len(na),len(cl),'FORMULA',int(np.floor(0.15*len(final_waters)/55.5+0.5)))
assert (len(final_waters),len(na),len(cl))==(3186,9,14)
# Flexible water preserves O-H bonds for a complete Amber parm7 topology.
sys=ff.createSystem(mod.topology,nonbondedMethod=PME,nonbondedCutoff=8.5*unit.angstrom,rigidWater=False)
boxvec=mod.topology.getPeriodicBoxVectors(); box6=[float(boxvec[0][0].value_in_unit(unit.angstrom)),float(boxvec[1][1].value_in_unit(unit.angstrom)),float(boxvec[2][2].value_in_unit(unit.angstrom)),90.,90.,90.]
structure=load_topology(mod.topology,system=sys,xyz=mod.positions,box=box6)
print('NONE_BONDS',sum(b.type is None for b in structure.bonds))
structure.save('/task/output/system.prmtop',overwrite=True)
structure.save('/task/output/system.inpcrd',format='rst7',overwrite=True)
pickle.dump((mod.topology,mod.positions,sys,box6),open('/task/output/final_built_state.pkl','wb'))
print('EXPORTED',len(structure.atoms),len(structure.residues),sum(a.charge for a in structure.atoms),box6)
