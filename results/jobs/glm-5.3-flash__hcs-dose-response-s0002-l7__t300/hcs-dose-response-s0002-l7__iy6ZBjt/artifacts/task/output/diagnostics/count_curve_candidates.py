import pandas as pd,numpy as np
from scipy.optimize import least_squares
c=pd.read_csv('/task/output/diagnostics/exact_counts.tsv',sep='\t')
c['row']=c.well.str[0].map(ord)-65;c['col']=c.well.str[1:].astype(int)-1
v=c[c.role=='vehicle'].reset_index(drop=True); y=v.n_cells.to_numpy(float)
def design(row,col,kind):
 row=np.asarray(row);col=np.asarray(col);rd=np.minimum(row,7-row);cd=np.minimum(col,11-col);x=(col-6)/5
 if kind=='rowcol2':
  X=np.zeros((len(row),8))
  for i,r in enumerate(row.astype(int)):X[i,r]=1
  return np.c_[X,x,x*x]
 if kind=='edge': return np.c_[np.ones_like(rd),rd,cd,rd*rd,cd*cd,rd*cd]
 if kind=='rowcollin':
  X=np.zeros((len(row),9))
  for i,r in enumerate(row.astype(int)):X[i,r]=1
  return np.c_[X,x]
models=['rowcol2','edge','rowcollin']
for kind in models:
 p=np.linalg.lstsq(design(v.row,v.col,kind),y,rcond=None)[0]
 print('\nMODEL',kind,'coef',p.round(2))
 for _,r in v.iterrows():
  pred=float(design(np.array([r.row]),np.array([r.col]),kind)@p)
  print(f'{r.well}:{r.n_cells}/{pred:.1f}={r.n_cells/pred:.2f}',end=' ')
 print()
 c[kind+'_pred']=design(c.row,c.col,kind)@p
 c[kind+'_resp']=c.n_cells/c[kind+'_pred']

for kind in models:
 print('\n### COUNT',kind)
 for cp,g in c[(c.role=='compound')&(~c.well.isin({'D06','D07','E07','E09','F07'}))].groupby('compound'):
  h=g.sort_values('dose',ascending=False);x=h.dose.to_numpy(float);yy=h[kind+'_resp'].to_numpy(float)
  best=None
  for b0 in [0,.25,.5,.75]:
   for h0 in [.5,1,2,3,5]:
    try:
     fun=lambda p:p[0]+(1-p[0])/(1+(x/max(p[1],1e-9))**p[2])-yy
     rr=least_squares(fun,[b0,max(np.min(x),1e-9),h0],bounds=([0,1e-9,.05],[.999,1000,100]),max_nfev=5000)
     if best is None or rr.cost<best[0]:best=(rr.cost,rr.x.copy())
    except:pass
  if best:
   sse,(b,e,hh)=best
   r=1/(1-2*b) if b<.5 else np.nan
   pot=e*r**(1/hh) if np.isfinite(r) and r>0 else np.nan
   print(f'{cp:10s} b={b:.2f} e={e:.3g} h={hh:.2f} GI50={pot:.3g} sse={sse:.3f} responses='+' '.join(f'{d:g}:{a:.2f}' for d,a in zip(x,yy)))
c.to_csv('/task/output/diagnostics/count_responses.tsv',sep='\t',index=False)
