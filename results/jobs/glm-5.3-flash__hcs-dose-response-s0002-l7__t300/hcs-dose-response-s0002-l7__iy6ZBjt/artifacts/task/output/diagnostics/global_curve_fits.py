import pandas as pd, numpy as np
from scipy.optimize import differential_evolution, least_squares
w=pd.read_csv('/task/output/well_measurements.tsv',sep='\t')
l=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')
d=w.merge(l[['well','compound','dose_uM','role']],on='well')
d=d[(d.qc_pass==1)&(d.role=='compound')]

def pred(x,p):
    b,le,h=p
    e=np.exp(le)
    z=np.clip(h*(np.log(np.maximum(x,1e-300))-le),-50,50)
    return b+(1-b)/(1+np.exp(z))

def global_fit(x,y):
    def sse(p): return float(np.sum((pred(x,p)-y)**2))
    r=differential_evolution(sse,[(0,.999999),np.log([1e-8,100]),(.05,100)],tol=1e-8,popsize=20,maxiter=500,polish=True,seed=123)
    b,le,h=r.x
    e=float(np.exp(le))
    pot=e*(1/(1-2*b))**(1/h) if b<.5 else np.nan
    return (float(b),e,float(h),float(pot),float(r.fun))
for cp,g in d.groupby('compound'):
    print('\n'+cp)
    for rd in ['response_n_cells','response_marker_per_cell']:
        g2=g.sort_values('dose_uM')
        x=g2.dose_uM.to_numpy(float); y=g2[rd].to_numpy(float)
        b,e,h,pot,sse=global_fit(x,y)
        print(f'{rd.replace("response_",""):6s} b={b:.4f} EC50={e:.6g} h={h:.3f} potency={pot:.6g} SSE={sse:.5f} responses='+' '.join(f'{a:g}:{z:.3f}' for a,z in zip(x,y)))
