import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import RBFInterpolator
m=pd.read_csv('/task/output/diagnostics/marker_candidates.tsv',sep='\t')
m['row']=m.well.str[0].map(ord)-65; m['col']=m.well.str[1:].astype(int)-1
v=m[m.role=='vehicle'].copy()

def poly_features(row,col,degree):
 out=[]
 for i in range(degree+1):
  for j in range(degree+1-i): out.append((row/7)**i*(col/11)**j)
 return np.column_stack(out)
# Leave-one-out CV on raw and log counts for polynomial and RBF models
for degree in [1,2,3,4,5]:
 errs=[]
 for i in v.index:
  tr=v.drop(i); te=v.loc[i]
  X=poly_features(tr.row.to_numpy(),tr.col.to_numpy(),degree)
  for target,use_log in [('n_cells',False),('n_cells',True)]:
   y=tr[target].to_numpy(float); yy=np.log(y) if use_log else y
   p=np.linalg.lstsq(X,yy,rcond=None)[0]
   pred=poly_features(np.array([te.row]),np.array([te.col]),degree)@p
   pred=np.exp(pred) if use_log else pred
   errs.append((target,use_log,pred[0]-te[target]))
 for t in ['n_cells']:
  e=[x[2] for x in errs if x[0]==t and not x[1]]; el=[x[2] for x in errs if x[0]==t and x[1]]
  print('poly',degree,'raw RMSE',np.sqrt(np.mean(np.array(e)**2)),'log RMSE ratio',np.sqrt(np.mean((np.array(el)/v.n_cells.to_numpy())**2)))
for degree in [1,2,3,4]:
 for eps in [0,1,5,10,20,50,100]:
  errs=[]
  for i in v.index:
   tr=v.drop(i); te=v.loc[i]
   rbf=RBFInterpolator(np.c_[tr.row/7,tr.col/11],tr.n_cells.to_numpy(float),kernel='thin_plate_spline',smoothing=eps,degree=degree)
   pred=rbf(np.array([[te.row/7,te.col/11]]))[0]
   errs.append((pred-te.n_cells)/te.n_cells)
  print('rbf',degree,'eps',eps,'relative RMSE',np.sqrt(np.mean(np.array(errs)**2)))
