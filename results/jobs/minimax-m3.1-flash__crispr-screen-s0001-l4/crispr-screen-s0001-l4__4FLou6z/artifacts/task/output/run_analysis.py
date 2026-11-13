#!/usr/bin/env python3
"""VAL-2 final count, QC, LFC, FDR-call, audit, and erratum analysis."""
import csv, gzip, json, os, subprocess
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import norm

INPUT='/task/input'; OUTPUT='/task/output'
SCAFFOLD='GTTTTAGAGCTA'
COMP=str.maketrans('ACGTN','TGCAN')
lib=pd.read_csv(f'{INPUT}/library.tsv',sep='\t',dtype=str)
samples=pd.read_csv(f'{INPUT}/samples.tsv',sep='\t',dtype=str)
reference=pd.read_csv(f'{INPUT}/reference_genes.tsv',sep='\t',dtype=str)
assert list(lib.columns)==['sgRNA','sequence','gene']
assert lib.sgRNA.is_unique and lib.sequence.is_unique and len(lib)==1000
seq_to_idx={s:i for i,s in enumerate(lib.sequence)}
counts=np.zeros((len(lib),len(samples)),dtype=np.int64)
qc_samples={}
for j,row in samples.iterrows():
    path=f'{INPUT}/fastq/{row.fastq_file}'
    assigned=Counter(); total=0
    with gzip.open(path,'rt') as fh:
        while True:
            h=fh.readline()
            if not h: break
            seq=fh.readline().strip(); fh.readline(); fh.readline(); total+=1
            # LIB15 is the NovaSeq Read 2 file. Reverse-complement it to sense orientation.
            if row.sample_id=='LIB15':
                seq=seq.translate(COMP)[::-1]
            # The spacer is immediately before the sgRNA scaffold. This single rule covers
            # all eight P5 stagger variants and tolerates sequencing errors in scaffold/constant.
            p=seq.find(SCAFFOLD)
            if p>=20:
                spacer=seq[p-20:p]
                k=seq_to_idx.get(spacer)
                if k is not None:
                    assigned[k]+=1
    counts[:,j]=[assigned[i] for i in range(len(lib))]
    qc_samples[row.sample_id]={'total_reads':int(total),'mapped_reads':int(counts[:,j].sum()),'arm':row.arm}

# Counts output (raw integer counts; one row per current-release library guide).
out=lib[['sgRNA','gene']].copy()
for j,sid in enumerate(samples.sample_id): out[sid]=counts[:,j]
out.to_csv(f'{OUTPUT}/counts.tsv',sep='\t',index=False)

# Profile-supported day-14 arms: LIB15 clusters with LIB11/LIB12. The known reference modifier
# signs orient the two clusters: vehicle=LIB11,LIB12,LIB15; drug=LIB10,LIB13,LIB14.
actual_arm={'LIB8':'T0','LIB9':'T0','LIB10':'drug','LIB11':'vehicle','LIB12':'vehicle','LIB13':'drug','LIB14':'drug','LIB15':'vehicle'}
for sid in qc_samples: qc_samples[sid]['arm']=actual_arm[sid]
with open(f'{OUTPUT}/qc.json','w') as f:
    json.dump({'samples':qc_samples},f,indent=2,sort_keys=False)

# Normalize depth using non-targeting median-ratio factors, then center each contrast on NTC=0.
ids=samples.sample_id.tolist(); gene=lib.gene.to_numpy(); ntc=gene=='NON_TARGETING'
C=counts.astype(float)
with np.errstate(all='ignore'):
    lm=np.log(np.where(C>0,C,np.nan))
    sf=np.exp(np.nanmedian(lm[ntc]-np.nanmean(lm[ntc],axis=1,keepdims=True),axis=0))
N=C/sf
t0=[ids.index(x) for x in ['LIB8','LIB9']]
veh=[ids.index(x) for x in ['LIB11','LIB12','LIB15']]
drug=[ids.index(x) for x in ['LIB10','LIB13','LIB14']]
def lfc(a,b):
    x=np.log2((N[:,a].mean(axis=1)+0.5)/(N[:,b].mean(axis=1)+0.5))
    return x-np.median(x[ntc])
fit=lfc(veh,t0); drg=lfc(drug,veh)

# Run MAGeCK for 5% FDR calls on the correctly grouped samples. Its gene-level neg/pos FDRs
# control the call set; required LFCs above remain the reported effect estimates.
control_file=f'{OUTPUT}/control_sgrna.txt'
with open(control_file,'w') as f:
    for sid in lib.loc[lib.gene=='NON_TARGETING','sgRNA']: f.write(sid+'\n')
prefix=f'{OUTPUT}/mageck_drug'
for p in [prefix+'.gene_summary.txt',prefix+'.sgrna_summary.txt',prefix+'.normalized_gene_summary.txt']:
    if os.path.exists(p): os.remove(p)
cmd=['mageck','test','-k',f'{OUTPUT}/counts.tsv','-t',','.join(ids[i] for i in drug),'-c',','.join(ids[i] for i in veh),'--control-sgrna',control_file,'-n',prefix]
subprocess.run(cmd,check=True,stdout=open(f'{OUTPUT}/mageck.log','w'),stderr=subprocess.STDOUT)
mageck=pd.read_csv(prefix+'.gene_summary.txt',sep='\t')

rows=[]
for g in sorted(set(gene[~ntc])):
    m=gene==g
    lf=float(np.median(fit[m])); ld=float(np.median(drg[m]))
    mr=mageck.loc[mageck.id==g]
    if len(mr)!=1: raise RuntimeError(f'MAGeCK missing gene {g}')
    mr=mr.iloc[0]
    neg_fdr=float(mr['neg|fdr']); pos_fdr=float(mr['pos|fdr'])
    if neg_fdr<=0.05 and (pos_fdr>0.05 or ld<0): cls='sensitiser'
    elif pos_fdr<=0.05 and (neg_fdr>0.05 or ld>0): cls='resistance'
    else: cls='none'
    rows.append((g,lf,ld,neg_fdr,pos_fdr,cls))
gs=pd.DataFrame(rows,columns=['gene','lfc_fitness','lfc_drug','fdr_sensitiser','fdr_resistance','drug_class'])
gs.to_csv(f'{OUTPUT}/gene_summary.tsv',sep='\t',index=False,float_format='%.6f')

# Release criterion validation.
ref=reference.merge(gs,on='gene',how='left')
mod_ok=int((ref.loc[ref.expected_class!='none','drug_class']==ref.loc[ref.expected_class!='none','expected_class']).sum())
none_fp=int((ref.loc[ref.expected_class=='none','drug_class']!='none').sum())
if mod_ok<9 or none_fp>1: raise RuntimeError(f'Release criterion failed: modifiers {mod_ok}/10, none false positives {none_fp}/30')

# Section 6 checklist, each statement judged independently.
audit=[
 ('library_release','defective'),
 ('sample_columns','ok'),
 ('stagger_reads','ok'),
 ('read2_library','defective'),
 ('arm_labels','defective'),
 ('control_normalisation','defective'),
 ('drug_baseline','ok'),
 ('contrast_sign','ok'),
]
pd.DataFrame(audit,columns=['item','verdict']).to_csv(f'{OUTPUT}/audit.tsv',sep='\t',index=False)

# Erratum for every mailed hit, in mailed order.
hits=pd.read_csv(f'{INPUT}/handover/hits_mailed.tsv',sep='\t')
classmap=gs.set_index('gene').drug_class.to_dict()
err=pd.DataFrame([(g,'confirmed' if classmap.get(g,'none') in ('sensitiser','resistance') else 'retracted') for g in hits.gene],columns=['gene','verdict'])
err.to_csv(f'{OUTPUT}/erratum.tsv',sep='\t',index=False)

print('release',mod_ok,'/10 modifiers;',none_fp,'/30 none false positives')
print(gs.drug_class.value_counts().to_dict())
print('mapped', {k:v['mapped_reads'] for k,v in qc_samples.items()})
