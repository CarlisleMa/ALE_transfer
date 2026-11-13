import pickle, numpy as np
from openmm import unit, Vec3
from openmm.app import Modeller, ForceField, PME
with open('/task/output/solvated_noions_state.pkl','rb') as f: top,pos=pickle.load(f)
ff=ForceField('amber14-all.xml')
sys=ff.createSystem(top,nonbondedMethod=PME,nonbondedCutoff=8.5*unit.angstrom)
charges=[]
for f_ in sys.getForces():
 if f_.getName()=='NonbondedForce':
  for i in range(sys.getNumParticles()):
   q,sigma,eps=f_.getParticleParameters(i); charges.append(q.value_in_unit(unit.elementary_charge))
protein_n=sum(1 for r in top.residues() if r.name!='HOH')
water_n=sum(1 for r in top.residues() if r.name=='HOH')
print('PROTEIN_RES',protein_n,'WATERS',water_n,'SOLUTE_CHARGE',sum(charges[:818]),'TOTAL',sum(charges))
print('SALT_ROUND',int(np.floor(0.15*water_n/55.5+0.5)))
# Ask OpenMM to place the exact desired ions in the same box on a throwaway copy.
# OpenMM's routine replaces waters; coordinates will be transplanted afterward.
tmp=Modeller(top,pos)
tmp._addIons(ff,water_n,{r:tmp.positions[next(a for a in r.atoms() if a.name=='O').index] for r in top.residues() if r.name=='HOH'},positiveIon='Na+',negativeIon='Cl-',ionicStrength=0*unit.molar,neutralize=True,residueTemplates={})
# Above neutralization only would give 5 Cl; then manually placement not exposed. Instead use ionic strength hack to place 14 Cl and 9 Na by setting custom charge? Print available ion residue coords.
ions=[r for r in tmp.topology.residues() if r.name in ('NA','CL','Na+','Cl-')]
print('TMP_IONS',[(r.name,r.id) for r in ions])
