from collections import defaultdict
src='/task/input/structure.pdb'; out='/task/output/protein_clean.pdb'
# Select protein ATOM records only. Select the highest-occupancy alternate per atom.
groups=defaultdict(list)
for line in open(src):
 if not line.startswith('ATOM  ') or line[21] != 'P': continue
 key=(line[22:26],line[12:16])
 groups[key].append(line)
selected=[]
for key,vals in groups.items():
 vals.sort(key=lambda x: (-float(x[54:60]), x[16]))
 selected.append(vals[0])
selected.sort(key=lambda x:(int(x[22:26]), x[26]))
# Build two chains at the unresolved gap, keeping deposited coordinates.
with open(out,'w') as f:
 f.write('REMARK  CLEANED RESOLVED PROTEIN; CHAIN A 18-36, CHAIN B 41-73\n')
 serial=1; chain='A'
 for line in selected:
  res=int(line[22:26]); ins=line[26]
  if res>=41 and chain=='A':
   f.write('TER\n'); chain='B'
  name=line[12:16].strip(); alt=line[16]; occ=line[54:60]; b=line[60:66]
  if res==67 and ins=='B':
   line=line[:17]+' '+line[17:20]+'MET'+line[21:]
   if name=='SE': line=line[:12]+' SD'+line[16:]
  outline=f'{line[:6]}{serial:5d} {line[11:16]} {chain}{line[21:26]}{ins}   {line[30:38]}{line[38:46]}{line[46:54]}{occ:>6}{b:>6}          {line[76:78]}\n'
  f.write(outline); serial+=1
 f.write('TER\nEND\n')
print(out,'atoms',serial-1)
