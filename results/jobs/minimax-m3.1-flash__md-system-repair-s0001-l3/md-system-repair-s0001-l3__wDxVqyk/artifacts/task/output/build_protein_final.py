from openmm.app import PDBFile, Topology, Modeller, ForceField, PDBFile as Writer
from openmm import unit, Vec3
from openmm.app.element import Element
import numpy as np, json
A=PDBFile('/task/output/fixed_A3.pdb'); B=PDBFile('/task/output/fixed_B3.pdb')
mapdata=json.load(open('/task/output/residue_mapping3.json'))
work_to_orig={x['work_num']:(str(x['old_num']),x['old_insertion']) for x in mapdata}

def arr(p): return np.array(p.value_in_unit(unit.angstrom),dtype=float)
def atom_pos(pdb,res_id,name):
 r=next(r for r in pdb.topology.residues() if str(r.id)==str(res_id))
 a=next(a for a in r.atoms() if a.name==name)
 return arr(pdb.positions[a.index])
def transform_points(ref_points,target_points):
    """Rigidly map three ref points onto three target points."""
    ref=np.asarray(ref_points,float); tar=np.asarray(target_points,float)
    r0=ref[0]; t0=tar[0]
    re1=ref[1]-r0; re2=ref[2]-r0
    te1=tar[1]-t0; te2=tar[2]-t0
    re1/=np.linalg.norm(re1); re2-=np.dot(re2,re1)*re1; re2/=np.linalg.norm(re2); re3=np.cross(re1,re2)
    te1/=np.linalg.norm(te1); te2-=np.dot(te2,te1)*te1; te2/=np.linalg.norm(te2); te3=np.cross(te1,te2)
    R=np.column_stack((te1,te2,te3))@np.column_stack((re1,re2,re3)).T
    return (np.asarray(ref_points,float)-r0)@R.T+t0
# Validated PDBFixer ACE/NME reference heavy-atom geometry.
ref_nme={'pre_CA':[4.853,4.614,0.0], 'pre_C':[4.713,6.129,0.0], 'pre_O':[3.601,6.653,0.0],
         'N':[5.846,6.835,0.0], 'CH3':[5.846,8.284,0.0]}
ref_ace={'post_N':[3.555,3.970,0.0], 'post_CA':[4.853,4.614,0.0], 'post_C':[4.713,6.129,0.0],
         'CH3':[2.000,2.090,0.0], 'C':[3.427,2.641,0.0], 'O':[4.391,1.877,0.0]}
tgt_nme=[atom_pos(A,'19','CA'),atom_pos(A,'19','C'),atom_pos(A,'19','O')]
tgt_ace=[atom_pos(B,'20','N'),atom_pos(B,'20','CA'),atom_pos(B,'20','C')]
nme_pts=transform_points([ref_nme['pre_CA'],ref_nme['pre_C'],ref_nme['pre_O']],tgt_nme)
ace_pts=transform_points([ref_ace['post_N'],ref_ace['post_CA'],ref_ace['post_C']],tgt_ace)
nme_all=transform_points([ref_nme[k] for k in ['pre_CA','pre_C','pre_O','N','CH3']],tgt_nme)
ace_all=transform_points([ref_ace[k] for k in ['post_N','post_CA','post_C','CH3','C','O']],tgt_ace)
# Values above correspond to transformed order: NME [CA,C,O,N,CH3], ACE [N,CA,C,CH3,C,O].
nme_N=nme_all[3]; nme_CH3=nme_all[4]
ace_CH3=ace_all[3]; ace_C=ace_all[4]; ace_O=ace_all[5]

top=Topology(); chain=top.addChain('P'); pos=[]; source_to_new={}; residue_by_source={}
def add_source_residue(pdb,oldres,skip=()):
    orig_id,ins=work_to_orig[int(oldres.id)]
    newres=top.addResidue(oldres.name,chain,orig_id,ins)
    residue_by_source[oldres]=newres
    for oldatom in oldres.atoms():
        if oldatom.name in skip: continue
        if oldres.name=='CYS' and oldatom.name=='HG': continue
        newatom=top.addAtom(oldatom.name,oldatom.element,newres)
        source_to_new[oldatom]=newatom
        pos.append(pdb.positions[oldatom.index])
    return newres
# Protein fragment A, then NME cap. Remove the temporary charged OXT from internal ILE36.
for r in A.topology.residues():
    skip={'OXT'} if str(r.id)=='19' else set()
    add_source_residue(A,r,skip)
# Map source bonds later, excluding atoms deliberately removed.
# NME cap.
nme=top.addResidue('NME',chain,'NME',' ')
a_nme_n=top.addAtom('N',Element.getBySymbol('N'),nme); pos.append(unit.Quantity(nme_N/10.0,unit.nanometer))
a_nme_ch3=top.addAtom('CH3',Element.getBySymbol('C'),nme); pos.append(unit.Quantity(nme_CH3/10.0,unit.nanometer))
# ACE cap, then protein fragment B. Remove temporary N-terminal H2/H3 from internal ASN41.
ace=top.addResidue('ACE',chain,'ACE',' ')
a_ace_ch3=top.addAtom('CH3',Element.getBySymbol('C'),ace); pos.append(unit.Quantity(ace_CH3/10.0,unit.nanometer))
a_ace_c=top.addAtom('C',Element.getBySymbol('C'),ace); pos.append(unit.Quantity(ace_C/10.0,unit.nanometer))
a_ace_o=top.addAtom('O',Element.getBySymbol('O'),ace); pos.append(unit.Quantity(ace_O/10.0,unit.nanometer))
for r in B.topology.residues():
    skip={'H2','H3'} if str(r.id)=='20' else set()
    add_source_residue(B,r,skip)
# Copy all source bonds that survived atom removal.
for pdb in (A,B):
    for a,b in pdb.topology.bonds():
        if a in source_to_new and b in source_to_new:
            top.addBond(source_to_new[a],source_to_new[b])
# Connect ACE to ASN41 and NME to ILE36.
def find_atom(res,name): return next(a for a in res.atoms() if a.name==name)
ile36=residue_by_source[next(r for r in A.topology.residues() if str(r.id)=='19')]
asn41=residue_by_source[next(r for r in B.topology.residues() if str(r.id)=='20')]
top.addBond(find_atom(ile36,'C'),a_nme_n)
top.addBond(a_nme_n,a_nme_ch3)
top.addBond(a_ace_ch3,a_ace_c)
top.addBond(a_ace_c,a_ace_o)
top.addBond(a_ace_c,find_atom(asn41,'N'))
# Add exactly the two cross-fragment disulfides; B already contains CYS47-CYS67A.
def source_res(pdb,rid): return next(r for r in pdb.topology.residues() if str(r.id)==str(rid))
for ai,bj in [('5','51'),('14','34')]:
    ra=residue_by_source[source_res(A,ai)]; rb=residue_by_source[source_res(B,bj)]
    top.addBond(find_atom(ra,'SG'),find_atom(rb,'SG'))
# Ensure no duplicate SG-SG bond was introduced.
sgbonds=[]
for a,b in top.bonds():
    if a.name=='SG' and b.name=='SG': sgbonds.append((a.residue.id,b.residue.id))
assert len(sgbonds)==3 and len(set(tuple(sorted(x)) for x in sgbonds))==3, sgbonds
# Rename bonded cysteines to CYX and remove thiol H.
for r in top.residues():
    if r.name=='CYS':
        sg=next((a for a in r.atoms() if a.name=='SG'),None)
        if sg is not None and any((x is sg and y.name=='SG') or (y is sg and x.name=='SG') for x,y in top.bonds()):
            r.name='CYX'
positions=unit.Quantity(np.asarray([[float(v) for v in q.value_in_unit(unit.angstrom)] for q in pos],dtype=float)/10.0,unit.nanometer)
mod=Modeller(top,positions)
ff=ForceField('amber14-all.xml')
# Existing protein hydrogens are retained; add the ACE/NME hydrogens.
# Add only the cap hydrogens explicitly; all protein hydrogens were already added by PDBFixer.
variants=[]
for r in mod.topology.residues():
    if r.name=='NME':
        variants.append([('H','N'),('HH31','CH3'),('HH32','CH3'),('HH33','CH3')])
    elif r.name=='ACE':
        variants.append([('HH31','CH3'),('HH32','CH3'),('HH33','CH3')])
    else:
        variants.append(None)
mod.addHydrogens(ff,variants=variants)
# Write prepared protein for inspection.
with open('/task/output/protein_final.pdb','w') as f: Writer.writeFile(mod.topology,mod.positions,f,keepIds=True)
# Save the exact topology/coordinates in a pickle for the next builder stage.
import pickle
with open('/task/output/protein_final_state.pkl','wb') as f:
    pickle.dump((mod.topology,mod.positions),f)
print('FINAL_PROTEIN',len(list(mod.topology.residues())),len(list(mod.topology.atoms())))
print('SG_BONDS',sgbonds)
print('CAPS',[(r.name,r.id,[a.name for a in r.atoms()]) for r in mod.topology.residues() if r.name in ('ACE','NME')])
print('TRUE_TERMI',[(r.name,r.id,[a.name for a in r.atoms()]) for r in [list(mod.topology.residues())[0],list(mod.topology.residues())[-1]]])
print('INTERNAL_ENDS',[(r.name,r.id,[a.name for a in r.atoms()]) for r in mod.topology.residues() if r.id in ('36','41','NME','ACE')])
