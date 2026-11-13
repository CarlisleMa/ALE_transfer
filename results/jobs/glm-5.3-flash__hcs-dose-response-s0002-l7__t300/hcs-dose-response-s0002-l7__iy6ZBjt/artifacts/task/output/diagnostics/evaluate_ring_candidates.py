import pandas as pd, numpy as np
m=pd.read_csv('/task/output/diagnostics/marker_candidates.tsv',sep='\t')
m['row']=m.well.str[0].map(ord)-65; m['col']=m.well.str[1:].astype(int)-1
qc_fail={'D06','D07','E07','E09','F07'}
choices=['ring3','ring4','ring5','ring6','ring8','mk_open15','mk_open25','mk_open40']
def design(row,col):
    row=np.asarray(row); col=np.asarray(col)
    X=np.zeros((len(row),8))
    for i,r in enumerate(row.astype(int)): X[i,r]=1
    x=(col-6)/5
    return np.c_[X,x,x*x]
for c in choices:
    v=m[(m.role=='vehicle')&m[c].notna()].reset_index(drop=True)
    y=np.log(np.maximum(v[c].to_numpy(float),1e-6))
    pred=np.empty(len(v))
    for i in range(len(v)):
        k=np.arange(len(v))!=i
        p=np.linalg.lstsq(design(v.row[k],v.col[k]),y[k],rcond=None)[0]
        pred[i]=float(design(np.array([v.row[i]]),np.array([v.col[i]]))@p)
    loo=np.sqrt(np.mean((pred-y)**2))
    p=np.linalg.lstsq(design(v.row,v.col),y,rcond=None)[0]
    m[c+'_pred']=np.exp(design(m.row,m.col)@p)
    m[c+'_resp']=m[c]/m[c+'_pred']
    print(f'\n### {c} LOO_LOG_RMSE={loo:.4f}')
    z=m[(m.role=='reference')&(~m.well.isin(qc_fail))]
    print('reference',' '.join(f'{r.well}:{r[c+"_resp"]:.2f}' for _,r in z.iterrows()))
    for cp,g in m[(m.role=='compound')&(~m.well.isin(qc_fail))].groupby('compound'):
        h=g.sort_values('dose',ascending=False)
        print(f'{cp:10s}',' '.join(f'{r.dose:g}:{r[c+"_resp"]:.2f}' for _,r in h.iterrows()))
