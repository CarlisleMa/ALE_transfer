import pickle, numpy as np
from openmm import unit
from openmm.app import Modeller, ForceField, PME
from parmed.openmm.topsystem import load_topology
FF_FILES=['amber14/protein.ff14SB.xml','amber14/tip3p.xml']
with open('/task/output/solvated_noions_state.pkl','rb') as f: top,pos=pickle.load(f)
mod=Modeller(top,pos); ff=ForceField(*FF_FILES)
# Complete oxygen-only solvent sites to TIP3P water; preserve protein hydrogens.
mod.addHydrogens(ff)
waters=[r for r in mod.topology.residues() if r.name=='HOH']
print('WATERS_COMPLETE',len(waters),'FIRST',[[a.name for a in r.atoms()] for r in waters[:2]])
# Determine the exact solute charge from the ff14SB nonbonded parameters.
probe=ff.createSystem(mod.topology,nonbondedMethod=PME,nonbondedCutoff=8.5*unit.angstrom)
nonbonded=next(f for f in probe.getForces() if f.getName()=='NonbondedForce')
protein_atoms=sum(len(list(r.atoms())) for r in mod.topology.residues() if r.name!='HOH')
charge=sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(protein_atoms))
print('PROTEIN_ATOMS',protein_atoms,'SOLUTE_CHARGE',charge)
# Place 5 Cl- for +5 solute and 9 NaCl pairs.  The internal routine uses 55.4 M;
# choose the concentration so its nearest-integer result is exactly 9 pairs.
n_waters=len(waters)
neutral=5
ionic=9*55.4/(n_waters-neutral)*unit.molar
water_pos={r:mod.positions[next(a for a in r.atoms() if a.name=='O').index] for r in waters}
mod._addIons(ff,n_waters,water_pos,positiveIon='Na+',negativeIon='Cl-',ionicStrength=ionic,neutralize=True,residueTemplates={})
final_waters=[r for r in mod.topology.residues() if r.name=='HOH']
na=[r for r in mod.topology.residues() if r.name=='NA']; cl=[r for r in mod.topology.residues() if r.name=='CL']
print('FINAL_COUNTS','water',len(final_waters),'Na',len(na),'Cl',len(cl),'salt_formula',int(np.floor(0.15*len(final_waters)/55.5+0.5)))
assert len(final_waters)==3186 and len(na)==9 and len(cl)==14
# Parameterize the complete system.
sys=ff.createSystem(mod.topology,nonbondedMethod=PME,nonbondedCutoff=8.5*unit.angstrom,rigidWater=True)
boxvec=mod.topology.getPeriodicBoxVectors()
box=[float(boxvec[0][0].value_in_unit(unit.angstrom)),float(boxvec[1][1].value_in_unit(unit.angstrom)),float(boxvec[2][2].value_in_unit(unit.angstrom)),90.0,90.0,90.0]
print('BOX',box,'SYSTEM_PARTICLES',sys.getNumParticles())
structure=load_topology(mod.topology,system=sys,xyz=mod.positions,box=box)
pickle.dump((mod.topology,mod.positions,sys),open('/task/output/final_built_state.pkl','wb'))
print('NONE_BONDS',[(i,b.atom1.name,b.atom1.residue.name,b.atom1.residue.number,b.atom2.name,b.atom2.residue.name,b.atom2.residue.number,b.type) for i,b in enumerate(structure.bonds) if b.type is None])
structure.save('/task/output/system.prmtop',overwrite=True)
structure.save('/task/output/system.inpcrd',format='rst7',overwrite=True)
pickle.dump((mod.topology,mod.positions),open('/task/output/final_built_state.pkl','wb'))
print('EXPORTED',len(structure.atoms),len(structure.residues),sum(a.charge for a in structure.atoms))
