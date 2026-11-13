from collections import defaultdict
src='/task/input/structure.pdb'
# Atom identity uses insertion code; altloc is not part of residue identity.
groups=defaultdict(list)
for idx,l in enumerate(open(src)):
 if l.startswith('ATOM  ') and l[21:22]=='P':
  groups[(l[22:26],l[26],l[12:16])].append((idx,l))
selected=[]
for key,vals in groups.items():
 selected.append(sorted(vals,key=lambda z:(-float(z[1][54:60]),z[1][16]))[0])
selected.sort()
# Correct residue identity and order: insertion code disambiguates 67, 67A, 67B.
res_order=[]
for _,l in selected:
 key=(int(l[22:26]),l[26])
 if key not in res_order: res_order.append(key)
work={key:i+1 for i,key in enumerate(res_order)}
print('RES ORDER',len(res_order),res_order)
for label,lo,hi in [('A',-999,36),('B',41,999)]:
 out=f'/task/output/fragment_{label}2.pdb'; serial=1; count=0
 with open(out,'w') as f:
  for _,l in selected:
   old=int(l[22:26]); ins=l[26]; key=(old,ins)
   if not (lo<=old<=hi): continue
   atom=l[12:16]; resn=l[17:20]
   if resn=='MSE':
    resn='MET'
    if atom.strip()=='SE': atom=' SD '
   line=f'{l[0:6]}{serial:5d} {atom:>4s}{l[16:17]}{resn:>3s} {label}{work[key]:4d}    {l[30:54]}{l[54:66]}\n'
   f.write(line); serial+=1
  f.write('TER\nEND\n')
 count=sum(1 for k in res_order if lo<=k[0]<=hi)
 print(label,out,count,'residues',serial-1,'atoms')
import json
json.dump([{'old_num':n,'old_insertion':i,'work_num':work[(n,i)]} for n,i in res_order],open('/task/output/residue_mapping2.json','w'),indent=2)
