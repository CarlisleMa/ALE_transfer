from collections import defaultdict
src='/task/input/structure.pdb'; out='/task/output/protein_clean3.pdb'
# Select protein atoms and highest-occupancy alternate conformer per atom name/residue.
all_atoms=[]; groups=defaultdict(list)
for idx,line in enumerate(open(src)):
    if line.startswith('ATOM  ') and line[21:22]=='P':
        all_atoms.append((idx,line))
        key=(line[22:26],line[26],line[12:16])
        groups[key].append((idx,line))
selected=[]
for key in sorted(groups):
    vals=sorted(groups[key],key=lambda z:(-float(z[1][54:60]),z[1][16]))
    selected.append(vals[0])
selected.sort()
# Build residue order directly from the deposited records. Sequence number + insertion code defines identity.
res_order=[]
for _,l in selected:
    key=(int(l[22:26]),l[26])
    if not res_order or res_order[-1]!=key: res_order.append(key)
# Map each deposited residue to a unique sequential OpenMM residue number.
resnum={key:i+1 for i,key in enumerate(res_order)}
# Split at unresolved 37-40.
with open(out,'w') as f:
    f.write('REMARK  RESOLVED PROTEIN ONLY; INTERNAL RESIDUES RENUMBERED FOR PARSING\n')
    serial=1; chain='A'; prev=None
    for _,l in selected:
        oldnum=int(l[22:26]); ins=l[26]; key=(oldnum,ins)
        if oldnum>=41 and chain=='A': f.write('TER\n'); chain='B'
        atom=l[12:16]; alt=l[16]; resn=l[17:20]
        if resn=='MSE':
            resn='MET'
            if atom.strip()=='SE': atom=' SD '
        # Unique sequential residue number and blank insertion code in working PDB.
        newline=f'{l[0:6]}{serial:5d} {atom:>4s}{alt}{resn:>3s} {chain}{resnum[key]:4d}    {l[30:54]}{l[54:66]}\n'
        f.write(newline); serial+=1
        prev=key
    f.write('TER\nEND\n')
# machine mapping
import json
mapping=[{'old_num':n,'old_insertion':i,'work_num':resnum[(n,i)]} for n,i in res_order]
json.dump(mapping,open('/task/output/residue_mapping.json','w'),indent=2)
print('residues',len(res_order),'atoms',len(selected),'mapping',mapping)
