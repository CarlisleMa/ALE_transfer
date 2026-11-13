import numpy as np, pandas as pd
m=pd.read_csv('/task/output/diagnostics/marker_candidates.tsv',sep='\t')
m['row']=m.well.str[0].map(ord)-65; m['col']=m.well.str[1:].astype(int)-1
v=m[m.role=='vehicle'].copy().reset_index(drop=True)
y=v.n_cells.to_numpy(float)

def poly(row,col,degree):
    return np.column_stack([(row/7)**i*(col/11)**j for i in range(degree+1) for j in range(degree+1-i)])

def rowcol(row,col,degree):
    # row fixed effects + centered global column polynomial
    D=np.zeros((len(row),8+(degree)))
    for k,r in enumerate(row): D[k,r]=1
    x=(col-6)/5
    for d in range(1,degree+1): D[:,8+d-1]=x**d
    return D

def edge(row,col):
    rd=np.minimum(row,7-row); cd=np.minimum(col,11-col)
    return np.c_[rd,cd,rd**2,cd**2,rd*cd]

models={f'poly{d}':lambda r,c,d=d:poly(r,c,d) for d in [1,2,3,4]}
models.update({f'rowcol{d}':lambda r,c,d=d:rowcol(r,c,d) for d in [1,2,3]})
models['edge']=lambda r,c:edge(r,c)
for name,fun in models.items():
    errs=[]
    for i in range(len(v)):
        keep=np.arange(len(v))!=i
        X=fun(v.row.to_numpy()[keep],v.col.to_numpy()[keep])
        p=np.linalg.lstsq(X,y[keep],rcond=None)[0]
        pred=fun(np.array([v.row[i]]),np.array([v.col[i]]))@p
        errs.append(pred[0]-y[i])
    rmse=np.sqrt(np.mean(np.square(errs)))
    X=fun(v.row.to_numpy(),v.col.to_numpy()); p=np.linalg.lstsq(X,y,rcond=None)[0]
    print(f'{name:10s} LOO_RMSE={rmse:5.1f} coefficients={np.round(p,1)}')
# Show global degree4 and row+quadratic predictions for all wells
for name in ['poly2','poly3','poly4','rowcol2']:
 fun=models[name]
 X=fun(v.row.to_numpy(),v.col.to_numpy()); p=np.linalg.lstsq(X,y,rcond=None)[0]
 pred=fun(m.row.to_numpy(),m.col.to_numpy())@p
 m[name]=pred
 print('\n',name)
 print(m.pivot_table(index='row',columns='col',values=name).round(1).to_string())
m.to_csv('/task/output/diagnostics/spatial_predictions.tsv',sep='\t',index=False)
