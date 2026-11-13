import numpy as np, pandas as pd
from scipy.optimize import least_squares
d=pd.read_csv('/task/output/explore_well_measurements.tsv',sep='\t')
d['row']=d.well.str[0]; d['col']=d.well.str[1:].astype(int)
qc=['C04','C05','E09','G06','G08']; u=d.query('well not in @qc').copy(); veh=u.query('role=="vehicle"').copy()

def positions(df):
    return np.c_[(df.row.map({r:i for i,r in enumerate('ABCDEFGH')})/7).values,(df.col-1)/11]
def baseline_count(metric,train,test):
    y=train[metric].values; g=y.mean(); rows=sorted(train.row.unique()); cols=sorted(train.col.unique())
    ri={r:i for i,r in enumerate(rows)}; ci={c:i for i,c in enumerate(cols)}
    re=np.zeros(len(rows)); ce=np.zeros(len(cols)); lam=0.5
    for _ in range(500):
        for i,r in enumerate(rows):
            ix=np.where(train.row.values==r)[0]
            re[i]=np.sum(y[ix]-g-ce[[ci[c] for c in train.col.values[ix]]])/(len(ix)+lam)
        for j,c in enumerate(cols):
            ix=np.where(train.col.values==c)[0]
            ce[j]=np.sum(y[ix]-g-re[[ri[r] for r in train.row.values[ix]]])/(len(ix)+lam)
    return g+np.array([re[ri[r]] for r in test.row])+np.array([ce[ci[c]] for c in test.col])
def baseline_marker(metric,train,test):
    from scipy.interpolate import RBFInterpolator
    f=RBFInterpolator(positions(train),train[metric].values,kernel='thin_plate_spline',smoothing=1.0)
    return np.maximum(f(positions(test)),np.finfo(float).eps)
def logistic(c,b,logc,h):
    c=np.asarray(c,float); ec=10**logc
    val=b+(1-b)/(1+(c/ec)**h)
    return np.where(c==0,1.0,val)
def fit(c,y):
    c=np.asarray(c,float); y=np.asarray(y,float); cmax=c[c>0].max()
    def fun(p): return logistic(c,*p)-y
    best=None
    for lb in np.linspace(0,.75,16):
        for logc in np.linspace(-6,2,17):
            for h in (.25,.5,1,2,4):
                z=least_squares(fun,[lb,logc,h],bounds=([0,-6,.05],[.999,6,100]),max_nfev=3000)
                if best is None or z.cost<best.cost: best=z
    b,logc,h=best.x; ec=10**logc
    cross=np.inf if b>=.5 else ec*(1/(1-2*b))**(1/h)
    return dict(bottom=b,ec50=ec,h=h,cross=cross,rss=2*best.cost,active=cross<=cmax)
for include_vehicle in (False,True):
    for metric,builder in [('n_area50',baseline_count),('r8_bgmode',baseline_marker),('r10_bgmode',baseline_marker)]:
        print(f'\n### include_vehicle={include_vehicle} metric={metric}')
        rows=[]
        for comp,g in u.query('role=="compound"').groupby('compound',sort=False):
            pred=builder(metric,veh,g)
            c=g.dose_uM.values; y=g[metric].values/pred
            if include_vehicle: c=np.r_[np.zeros(len(veh)),c]; y=np.r_[np.ones(len(veh)),y]
            z=fit(c,y)
            rows.append((comp,float(y.min()),float(y.max()),z['bottom'],z['h'],z['ec50'],z['cross'],'active' if z['active'] else 'inactive'))
        print(pd.DataFrame(rows,columns=['compound','minresp','maxresp','bottom','h','potency','halfvehicle','decision']).round(4).to_string(index=False))
