import pandas as pd, numpy as np
from scipy.optimize import least_squares
c=pd.read_csv('/task/output/diagnostics/exact_counts.tsv',sep='\t')
m=pd.read_csv('/task/output/diagnostics/fast_marker_candidates.tsv',sep='\t')
d=c.merge(m[['well','remote_2-8']],on='well')
d['row']=d.well.str[0].map(ord)-65; d['col']=d.well.str[1:].astype(int)-1
qc={'D06','D07','E07','E09','F07'}
v=d[d.role=='vehicle'].copy()
def design(row,col):
    row=np.asarray(row); col=np.asarray(col)
    X=np.zeros((len(row),8))
    for i,r in enumerate(row.astype(int)): X[i,r]=1
    z=(col-6)/5
    return np.c_[X,z,z*z]
for target in ['n_cells_exact','remote_2-8']:
    y=v[target].to_numpy(float)
    link=np.log(np.maximum(y,1e-6)) if target=='remote_2-8' else y
    p=np.linalg.lstsq(design(v.row,v.col),link,rcond=None)[0]
    pred=design(d.row,d.col)@p
    if target=='remote_2-8': pred=np.exp(pred)
    d['resp']=d[target]/pred
    print('\n###',target)
    for loss,scale in [('linear',None),('soft_l1',0.1),('soft_l1',0.2),('huber',0.1),('cauchy',0.1)]:
        print('LOSS',loss,scale)
        for cp,g in d[(d.role=='compound')&(~d.well.isin(qc))].groupby('compound'):
            h=g.sort_values('dose',ascending=False); x=h.dose.to_numpy(float); yy=h.resp.to_numpy(float)
            best=None
            for b0 in [0,.25,.5,.75]:
                for h0 in [.5,1,2,3,5]:
                    f=lambda q:q[0]+(1-q[0])/(1+(x/np.maximum(q[1],1e-9))**q[2])-yy
                    try:
                        kw={} if loss=='linear' else {'loss':loss,'f_scale':scale}
                        r=least_squares(f,[b0,max(np.min(x),1e-9),h0],bounds=([0,1e-9,.05],[.999,1000,100]),max_nfev=5000,**kw)
                        score=np.sum(r.fun**2)
                        if best is None or score<best[0]: best=(score,r.x.copy())
                    except Exception: pass
            if best:
                sse,(b,e,hh)=best
                ratio=2/(1-b)-1
                pot=e*ratio**(1/hh) if ratio>0 else np.nan
                active=np.isfinite(pot) and pot<=30
                print(f'{cp:10s} b={b:.2f} ec={e:.3g} h={hh:.2f} pot={pot:.3g} active={active} responses='+' '.join(f'{a:g}:{z:.2f}' for a,z in zip(x,yy)))
