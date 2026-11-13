import pandas as pd, numpy as np
from scipy.optimize import least_squares
m=pd.read_csv('/task/output/diagnostics/fast_marker_candidates.tsv',sep='\t')
m['row']=m.well.str[0].map(ord)-65; m['col']=m.well.str[1:].astype(int)-1
qc_fail={'D06','D07','E07','E09','F07'}
vals=[c for c in m if c.startswith(('raw_','clip_','remote_'))]

def design(row,col,degree=3):
    row=np.asarray(row); col=np.asarray(col)
    X=np.zeros((len(row),9))
    for i,r in enumerate(row.astype(int)): X[i,r]=1
    x=(col-6)/5
    return np.c_[X,x,x**2,x**3]
# LOO on vehicle for each marker candidate
errs={}
for c in vals:
    v=m[(m.role=='vehicle')&m[c].notna()].reset_index(drop=True)
    y=np.log(np.maximum(v[c].to_numpy(float),1e-6))
    pred=np.empty(len(v))
    for i in range(len(v)):
        k=np.arange(len(v))!=i
        p=np.linalg.lstsq(design(v.row[k],v.col[k]),y[k],rcond=None)[0]
        pred[i]=design(np.array([v.row[i]]),np.array([v.col[i]]))@p
    errs[c]=np.sqrt(np.mean((pred-y)**2))
print('LOO LOG RMSE')
for c,e in sorted(errs.items(),key=lambda z:z[1]): print(f'{e:.4f} {c}')
# Fit surface to all vehicle and normalize all wells
for c in vals:
    v=m[(m.role=='vehicle')&m[c].notna()]
    p=np.linalg.lstsq(design(v.row,v.col),np.log(np.maximum(v[c].to_numpy(float),1e-6)),rcond=None)[0]
    m[c+'_pred']=np.exp(design(m.row,m.col)@p)
    m[c+'_resp']=m[c]/m[c+'_pred']
print('\nNORMALIZED CONTROLS')
for c in vals:
 z=m[m.role=='reference'][['well',c+'_resp']]
 print(c, ' '.join(f'{r.well}:{r[c+"_resp"]:.2f}' for r in z.iterrows()))
print('\nCOMPOUND NORMALIZED DOSE SERIES')
for c in vals:
 print('\n###',c)
 for cp,g in m[(m.role=='compound')&(~m.well.isin(qc_fail))].groupby('compound'):
  h=g.sort_values('dose',ascending=False)
  print(f'{cp:10s}',' '.join(f'{r.dose:g}:{r[c+"_resp"]:.2f}' for r in h.iterrows()))
m.to_csv('/task/output/diagnostics/marker_candidates_normalized.tsv',sep='\t',index=False)
