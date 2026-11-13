from collections import defaultdict
src='/task/input/structure.pdb'; out='/task/output/protein_clean2.pdb'
groups=defaultdict(list)
idx=0
for line in open(src):
    if line.startswith('ATOM  ') and line[21:22]=='P':
        # Highest-occupancy alternate for each atom; deterministic A before B on ties.
        groups[(line[22:26],line[12:16])].append((idx,line)); idx+=1
selected=[]
for key in sorted(groups):
    vals=sorted(groups[key],key=lambda x:(-float(x[1][54:60]),x[1][16]))
    selected.append(vals[0])
# Keep deposited atom order from the input, not lexical order.
selected.sort(key=lambda x:x[0])
with open(out,'w') as f:
    f.write('REMARK  RESOLVED PROTEIN ONLY; UNRESOLVED 37-40 EXCLUDED\n')
    serial=1; chain='A'
    for _,line in selected:
        res=int(line[22:26]); ins=line[26]
        if res>=41 and chain=='A':
            f.write('TER\n'); chain='B'
        # Copy original fixed columns, changing only chain, residue name for MSE, and atom name SE->SD.
        atom=line[12:16]
        resn=line[17:20]
        if resn=='MSE': resn='MET'
        if resn=='MSE' and atom.strip()=='SE': atom=' SD '
        newline=line[:12]+atom+line[16:17]+resn+line[20:21]+chain+line[22:26]+ins+line[27:30]
        # Use original coordinate/occupancy fields; atom serial/line formatting corrected.
        newline=f'{line[0:6]}{serial:5d} {atom:>4s}{line[16:17]}{resn:>3s} {chain}{line[22:26]}{ins}   {line[30:54]}{line[54:66]}\n'
        f.write(newline); serial+=1
    f.write('TER\nEND\n')
print('WROTE',out,'atoms',serial-1)
