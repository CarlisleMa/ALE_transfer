import pandas as pd, numpy as np
from scipy.optimize import least_squares
w=pd.read_csv('/task/output/well_measurements.tsv',sep='\t')
l=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')
d=w.merge(l[['well','compound','dose_uM','role']],on='well')
d=d[(d.qc_pass==1)&(d.role=='compound')]
def model(x,p):
    b,e,h=p
    z=np.clip(h*np.log(np.maximum(x,1e-300)/max(e,1e-300)),-50,50)
    return b+(1-b)/(1+np.exp(z))
def fit(x,y,loss='linear',constrained=False):
    best=None
    xmax=x.max()
    bounds=([0,1e-12,.05],[.999999,1e6,100])
    starts=[]
    for b0 in [0,.25,.5,.75]:
        for h0 in [.5,1,2,3,5,10]:
            starts.append([b0,x.max()/2,h0])
            if not constrained: starts.append([b0,x.min()/2,h0])
    for p0 in starts:
        try:
            r=least_squares(lambda p:model(x,p)-y,p0,bounds=bounds,loss=loss,f_scale=.2,max_nfev=10000)
            score=float(np.sum(r.fun**2))
            if best is None or score<best[0]: best=(score,r.x.copy())
        except: pass
    b,e,h=best[1]
    pot=e*(1/(1-2*b))**(1/h) if b<.5 else np.nan
    return best[1],pot,best[0]
for cp,g in d.groupby('compound'):
    x=g.dose_uM.to_numpy(); print('\n',cp)
    for rd in ['response_marker_per_cell','response_n_cells']:
        y=g[rd].to_numpy()
        for loss in ['linear','soft_l1']:
            for constr in [False,True]:
                p,pot,sse=fit(x,y,loss,constr)
                print(f'{rd.replace("response_",""):6s} {loss:7s} constrained={str(constr):5s} b={p[0]:.3f} EC={p[1]:.4g} h={p[2]:.2f} pot={pot:.4g} SSE={sse:.3f}')
