from collections import defaultdict
src='/task/input/structure.pdb'
groups=defaultdict(list); order=[]
for idx,l in enumerate(open(src)):
 if l.startswith('ATOM  ') and l[21:22]=='P':
  groups[(l[22:26],l[26],l[12:16])].append((idx,l)); order.append((idx,l))
selected=[]
for key,vals in groups.items():
 selected.append(sorted(vals,key=lambda z:(-float(z[1][54:60]),z[1][16]))[0])
selected.sort()
# derive residue order and unique work numbering globally
res_order=[]
for _,l in selected:
 key=(int(l[22:26]),l[26])
 if not res_order or res_order[-1]!=key: res_order.append(key)
work={key:i+1 for i,key in enumerate(res_order)}
for label,lo,hi in [('A',-999,36),('B',41,999)]:
 out=f'/task/output/fragment_{label}.pdb'; serial=1
 with open(out,'w') as f:
  for _,l in selected:
   old=int(l[22:26]); ins=l[26]; key=(old,ins)
   if not (lo<=old<=hi): continue
   atom=l[12:16]; resn=l[17:20]
   if resn=='MSE':
    resn='MET'
    if atom.strip()=='SE': atom=' SD '
   # chain IDs unique in each fragment; work IDs globally unique.
   line=f'{l[0:6]}{serial:5d} {atom:>4s}{l[16:17]}{resn:>3s} {label}{work[key]:4d}    {l[30:54]}{l[54:66]}\n'
   f.write(line); serial+=1
  f.write('TER\nEND\n')
 print(label,out,sum(1 for k in res_order if lo<=k[0]<=hi),'residues')
