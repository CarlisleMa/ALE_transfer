import numpy as np, pandas as pd
from scipy.optimize import least_squares
c=pd.read_csv('/task/output/diagnostics/exact_counts.tsv',sep='\t')
m=pd.read_csv('/task/output/diagnostics/fast_marker_candidates.tsv',sep='\t')
d=c.merge(m[['well','remote_2-6','remote_2-8','g20_2-6']],on='well')
d['row']=d.well.str[0].map(ord)-65; d['col']=d.well.str[1:].astype(int)-1
qc={'D06','D07','E07','E09','F07'}

def feat(row,col,kind):
    row=np.asarray(row,float); col=np.asarray(col,float)
    r=row/7; z=col/11; rd=np.minimum(row,7-row)/3.5; cd=np.minimum(col,11-col)/5.5
    if kind=='lin': return np.c_[np.ones_like(r),r,z]
    if kind=='quad': return np.c_[np.ones_like(r),r,z,r*r,z*z,r*z]
    if kind=='cubic':
        return np.column_stack([r**i*z**j for i in range(4) for j in range(4-i)])
    if kind=='rowquad':
        X=np.zeros((len(row),8))
        for i,x in enumerate(row.astype(int)): X[i,x]=1
        return np.c_[X,(z-.5),(z-.5)**2]
    if kind=='rowcubic':
        X=np.zeros((len(row),8))
        for i,x in enumerate(row.astype(int)): X[i,x]=1
        zz=(z-.5)
        return np.c_[X,zz,zz**2,zz**3]
    if kind=='edge': return np.c_[np.ones_like(rd),rd,cd,rd*rd,cd*cd,rd*cd]

def loess_features(qrow,qcol,trow,tcol,deg=1,band=2.5):
    # normalized plate coordinates; tricube weighted local polynomial
    q=np.c_[np.asarray(qrow,float)/7,np.asarray(qcol,float)/11]
    t=np.c_[np.asarray(trow,float)/7,np.asarray(tcol,float)/11]
    out=[]
    for p in q:
        dist=np.sqrt(np.sum((t-p)**2,axis=1)); h=max(band,np.partition(dist,1)[1]*1.15)
        w=np.clip((1-(dist/h)**3)**3,0,None); w=w/w.sum()
        # weighted centered monomials
        rr=t[:,0]-p[0]; cc=t[:,1]-p[1]
        if deg==1: M=np.c_[np.ones_like(rr),rr,cc]
        else: M=np.c_[np.ones_like(rr),rr,cc,rr*rr,cc*cc,rr*cc]
        root=np.sqrt(w)[:,None]
        out.append((M,root,p))
    return out

def fit_predict(kind,trr,trc,ter,qrr,qrc):
    if kind.startswith('loess'):
        deg=1 if kind=='loess1' else 2
        fits=loess_features(qrr,qrc,trr,trc,deg)
        yr=np.asarray(ter,float)
        pred=[]
        for M,root,p in fits:
            beta=np.linalg.lstsq(M*root,yr[:,None]*root,rcond=None)[0]
            featv=np.r_[1.,0.,0.] if deg==1 else np.r_[1.,0.,0.,0.,0.,0.]
            pred.append(float(featv@beta))
        return np.asarray(pred)
    X=feat(trr,trc,kind); y=np.asarray(ter,float)
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    return feat(qrr,qrc,kind)@beta

models=['lin','quad','cubic','rowquad','rowcubic','edge','loess1','loess2']
for target in ['n_cells_exact','remote_2-6','remote_2-8','g20_2-6']:
    v=d[(d.role=='vehicle')&d[target].notna()].reset_index(drop=True)
    y=v[target].to_numpy(float)
    use=np.log(np.maximum(y,1e-6)) if target!='n_cells_exact' else y
    print('\nTARGET',target)
    for kind in models:
        pred=fit_predict(kind,v.row,v.col,use,v.row,v.col)
        rmse=np.sqrt(np.mean((pred-use)**2))
        allpred=fit_predict(kind,v.row,v.col,use,d.row,d.col)
        if target!='n_cells_exact': allpred=np.exp(allpred)
        key='pred_'+kind
        d[key]=allpred
        d['resp_'+target+'_'+kind]=d[target]/d[key]
        print(f'{kind:9s} LOO_RMSE={rmse:.4f}',end='  ')
    print()
# fit curves for each model/readout

def fit_curve(x,y):
    best=None
    for b0 in [0,.25,.5,.75]:
        for h0 in [.5,1,2,3,5]:
            try:
                f=lambda p:p[0]+(1-p[0])/(1+(x/np.maximum(p[1],1e-9))**p[2])-y
                r=least_squares(f,[b0,max(np.min(x),1e-9),h0],bounds=([0,1e-9,.05],[.999,1000,100]),max_nfev=5000)
                if best is None or r.cost<best[0]: best=(r.cost,r.x.copy())
            except Exception: pass
    if best is None:return None
    sse,(b,e,h)=best
    rr=2/(1-b)-1
    pot=e*rr**(1/h) if rr>0 else np.nan
    return b,e,h,pot,sse
for target in ['n_cells_exact','remote_2-6','remote_2-8','g20_2-6']:
    print('\nCURVES',target)
    for kind in models:
        print(' MODEL',kind)
        for cp,g in d[(d.role=='compound')&(~d.well.isin(qc))].groupby('compound'):
            h=g.sort_values('dose',ascending=False)
            x=h.dose.to_numpy(float); y=h['resp_'+target+'_'+kind].to_numpy(float)
            z=fit_curve(x,y)
            if z:
                b,e,hh,pot,sse=z
                print(f'{cp:10s} b={b:.2f} ec={e:.3g} h={hh:.2f} pot={pot:.3g} sse={sse:.3f}')
