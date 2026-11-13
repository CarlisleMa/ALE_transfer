import pandas as pd, numpy as np
from scipy.optimize import least_squares
c=pd.read_csv('/task/output/diagnostics/exact_counts.tsv',sep='\t')
m=pd.read_csv('/task/output/diagnostics/polynomial_marker.tsv',sep='\t')
d=c.merge(m[['well','deg2_2-8']],on='well')
d['row']=d.well.str[0].map(ord)-65; d['col']=d.well.str[1:].astype(int)-1
qc={'D06','D07','E07','E09','F07'}

def piecewise_extrap(v,col):
    v=v.sort_values('col')
    x=v.col.to_numpy(float); y=v.target.to_numpy(float)
    if len(x)==1:return np.full(len(col),y[0])
    out=np.interp(col,x,y)
    lo=col<x[0]; hi=col>x[-1]
    if lo.any():
        slope=(y[1]-y[0])/(x[1]-x[0]); out[lo]=y[0]+slope*(col[lo]-x[0])
    if hi.any():
        slope=(y[-1]-y[-2])/(x[-1]-x[-2]); out[hi]=y[-1]+slope*(col[hi]-x[-1])
    return out

def rowquad(row,col):
    row=np.asarray(row,int); z=(np.asarray(col)-6)/5
    X=np.zeros((len(row),8))
    for i,r in enumerate(row): X[i,r]=1
    return np.c_[X,z,z*z]

def fitpot(x,y):
    best=None; x0=float(np.exp(np.mean(np.log(x))))
    def model(p):
        arg=np.clip(p[2]*np.log(np.maximum(x,1e-300)/max(p[1],1e-300)),-50,50)
        return p[0]+(1-p[0])/(1+np.exp(arg))-y
    for b0 in [0,.25,.5,.75]:
        for h0 in [.5,1,2,3,5,10]:
            try:
                r=least_squares(model,[b0,x0,h0],bounds=([0,1e-12,.05],[.999999,1e6,100]),max_nfev=10000)
                s=float(np.sum(r.fun**2))
                if best is None or s<best[0]:best=(s,r.x.copy())
            except: pass
    b,e,h=best[1]; pot=e*(1/(1-2*b))**(1/h) if b<.5 else np.nan
    return b,e,h,pot,best[0]
for target in ['n_cells_exact','deg2_2-8']:
    print('\n###',target)
    # global rowquad
    v=d[(d.role=='vehicle')&d[target].notna()]
    y=v[target].to_numpy(float)
    use=np.log(np.maximum(y,1e-12)) if target!='n_cells_exact' else y
    p=np.linalg.lstsq(rowquad(v.row,v.col),use,rcond=None)[0]
    pred=rowquad(d.row,d.col)@p
    if target!='n_cells_exact': pred=np.exp(pred)
    d['global_pred']=pred
    # row-specific interpolation
    preds=[]
    for _,r in d.iterrows():
        vv=v[v.row==r.row]
        preds.append(piecewise_extrap(vv.assign(target=vv[target]),np.array([r.col]))[0])
    d['interp_pred']=preds
    for method in ['global_pred','interp_pred']:
        d['resp']=d[target]/d[method]
        print('\nMETHOD',method)
        for cp,g in d[(d.role=='compound')&(~d.well.isin(qc))].groupby('compound'):
            h=g.sort_values('dose'); x=h.dose.to_numpy(); yy=h.resp.to_numpy()
            b,e,hh,pot,sse=fitpot(x,yy)
            inrange=np.isfinite(pot) and x.min()<=pot<=x.max()
            print(f'{cp:10s} pot={pot:.5g} inrange={inrange} responses='+' '.join(f'{a:g}:{z:.2f}' for a,z in zip(x,yy)))
