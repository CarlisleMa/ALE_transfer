import numpy as np,pandas as pd
from scipy.interpolate import RBFInterpolator
d=pd.read_csv('/task/output/explore_well_measurements.tsv',sep='\t')
d['row']=d.well.str[0]; d['col']=d.well.str[1:].astype(int)
qc=['C04','C05','E09','G06','G08']; u=d.query('well not in @qc').copy(); veh=u.query('role=="vehicle"').copy()
# positions normalized independently for interpolation
pos=lambda df:np.c_[(df.row.map({r:i for i,r in enumerate('ABCDEFGH')})/7),(df.col-1)/11]
def predictions(metric,train,test,kind):
 y=train[metric].values
 if kind=='global': return np.repeat(np.mean(y),len(test))
 if kind=='median': return np.repeat(np.median(y),len(test))
 if kind=='row':
  mp=train.groupby('row')[metric].mean(); return test.row.map(mp).fillna(train[metric].mean()).values
 if kind=='col':
  mp=train.groupby('col')[metric].mean(); return test.col.map(mp).fillna(train[metric].mean()).values
 if kind=='rowcat_collinear':
  # intercept + row indicators + linear column
  rows=sorted(train.row.unique()); X=np.c_[np.ones(len(train)),np.array([[r==q for q in rows] for r in train.row],float),train.col]
  beta=np.linalg.lstsq(X,y,rcond=None)[0]
  Xt=np.c_[np.ones(len(test)),np.array([[r==q for q in rows] for r in test.row],float),test.col]
  return Xt@beta
 if kind=='additive_cat':
  # alternating least squares categorical row+column with ridge regularization
  rows=sorted(train.row.unique()); cols=sorted(train.col.unique())
  ri={r:i for i,r in enumerate(rows)}; ci={c:i for i,c in enumerate(cols)}
  g=np.mean(y); re=np.zeros(len(rows)); ce=np.zeros(len(cols)); lam=1.0
  for _ in range(100):
   for i,r in enumerate(rows):
    ix=np.where(train.row.values==r)[0]
    re[i]=np.sum(y[ix]-g-ce[([ci[c] for c in train.col.values[ix]])])/(len(ix)+lam)
   for j,c in enumerate(cols):
    ix=np.where(train.col.values==c)[0]
    ce[j]=np.sum(y[ix]-g-re[([ri[r] for r in train.row.values[ix]])])/(len(ix)+lam)
  return g+np.array([re[ri[r]] for r in test.row])+np.array([ce[ci[c]] for c in test.col])
 if kind=='medianpolish':
  tab=train.pivot(index='row',columns='col',values=metric)
  g=tab.stack().median(); re=pd.Series(0.,index=tab.index); ce=pd.Series(0.,index=tab.columns)
  for _ in range(200):
   re += (tab.sub(g+ce,axis=1)).median(axis=1)-re
   g=(tab.sub(re,axis=0).sub(ce)).stack().median()
   ce += (tab.sub(re,axis=0)).median(axis=0)-ce
  return np.array([g+re.loc[r]+ce.loc[c] for r,c in zip(test.row,test.col)])
 if kind.startswith('rbf'):
  smoothing=float(kind.split('_')[1]); f=RBFInterpolator(pos(train),y,kernel='thin_plate_spline',smoothing=smoothing)
  return f(pos(test))
for metric in ['n_area40','r8_bgmode','r10_bgmode']:
 print('\n###',metric)
 kinds=['global','median','row','col','rowcat_collinear','additive_cat','medianpolish','rbf_0','rbf_1','rbf_10','rbf_100']
 for kind in kinds:
  errs=[]
  # If a categorical level is missing when left out, model may still use other levels; catch.
  for i in veh.index:
   tr=veh.drop(index=i); te=veh.loc[[i]]
   try:
    pred=predictions(metric,tr,te,kind); errs.extend((te[metric].values-pred).tolist())
   except Exception as e: errs.append(np.nan)
  rmse=np.sqrt(np.nanmean(np.array(errs)**2))
  pred=predictions(metric,veh,veh,kind)
  print(f'{kind:20s} LOO_RMSE={rmse:8.2f} residual_sd={np.std(veh[metric]-pred):8.2f}')
