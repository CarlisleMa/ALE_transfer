import numpy as np,pandas as pd
from scipy.optimize import differential_evolution, least_squares
d=pd.read_csv('/task/output/explore_well_measurements.tsv',sep='\t')
d['row']=d.well.str[0]; d['col']=d.well.str[1:].astype(int)
qc=['C04','C05','E09','G06','G08']; u=d.query('well not in @qc').copy(); veh=u.query('role=="vehicle"').copy()

def positions(df):
 return np.c_[(df.row.map({r:i for i,r in enumerate('ABCDEFGH')})/7).values,(df.col-1)/11]
def baseline_count(metric,train,test):
 y=train[metric].values; g=y.mean(); rows=sorted(train.row.unique()); cols=sorted(train.col.unique())
 ri={r:i for i,r in enumerate(rows)}; ci={c:i for i,c in enumerate(cols)}
 re=np.zeros(len(rows)); ce=np.zeros(len(cols)); lam=.5
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
 c=np.asarray(c,float)
 ec=np.exp(logc)
 x=np.full(c.shape,0.0) if False else c/ec
 val=b+(1-b)/(1+np.where(c==0,0,x)**h)
 return np.where(c==0,1.0,val)
def fit(c,y):
 c=np.asarray(c,float); y=np.asarray(y,float); cmax=c[c>0].max()
 loglo=np.log(1e-6); loghi=np.log(1e6)
 def fun(p): return logistic(c,*p)-y
 def objective(p):
  z=fun(p); return float(np.dot(z,z))
 de=differential_evolution(objective,[(0,.999),(loglo,loghi),(.05,100)],tol=1e-9,polish=True,seed=12345,maxiter=1000,popsize=25)
 z=least_squares(fun,de.x,bounds=([0,loglo,.05],[.999,loghi,100]),max_nfev=10000)
 b,logc,h=z.x; ec=np.exp(logc)
 if b>=.5: cross=np.inf
 else: cross=ec*(1-2*b)**(1/h)
 return dict(bottom=b,ec50=ec,h=h,cross=cross,rss=2*z.cost,active=cross<=cmax,cmax=cmax)
for include_vehicle in (False,True):
 for metric,builder in [('n_area40',baseline_count),('r8_bgmode',baseline_marker),('r10_bgmode',baseline_marker)]:
  print(f'\n### include_vehicle={include_vehicle} metric={metric} baseline={builder.__name__}')
  rows=[]
  for comp,g in u.query('role=="compound"').groupby('compound',sort=True):
   pred=builder(metric,veh,g)
   allc=np.r_[np.zeros(len(veh)) if include_vehicle else [],g.dose_uM.values]
   ally=np.r_[veh[metric].values/pred[0]*0+1 if include_vehicle else [],g[metric].values/pred]
   z=fit(allc,ally)
   resp=g[metric].values/pred
   rows.append((comp,min(resp),max(resp),z['bottom'],z['h'],z['ec50'],z['cross'],z['cmax'],'active' if z['active'] else 'inactive'))
  print(pd.DataFrame(rows,columns=['compound','minresp','maxresp','bottom','h','ec50','halfvehicle','maxtested','decision']).round(4).to_string(index=False))
