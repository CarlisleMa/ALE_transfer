import pandas as pd,numpy as np
from scipy.stats import ttest_1samp
D=pd.read_csv('/task/output/counts.tsv',sep='\t'); ids=D.columns[2:].tolist(); L=pd.read_csv('/task/input/library.tsv',sep='\t'); R=pd.read_csv('/task/input/reference_genes.tsv',sep='\t'); gene=L.gene.to_numpy(); C=D[ids].to_numpy(float); ntc=gene=='NON_TARGETING'; veh=[3,4,7]; drg=[2,5,6]
with np.errstate(all='ignore'):
 lm=np.log(np.where(C>0,C,np.nan)); sf=np.exp(np.nanmedian(lm[ntc]-np.nanmean(lm[ntc],axis=1,keepdims=True),axis=0))
N=C/sf
def lfc(a,b):
 x=np.log2((N[:,a].mean(1)+.5)/(N[:,b].mean(1)+.5)); return x-np.median(x[ntc])
fit=lfc(veh,[0,1]); d=lfc(drg,veh); mu=np.median(d[ntc]); sd=np.std(d[ntc],ddof=1)
rows=[]
for g in sorted(set(gene[~ntc])):
 m=gene==g; vals=d[m]; p=ttest_1samp(vals,0).pvalue; rows.append((g,np.median(fit[m]),np.median(vals),p,np.std(vals,ddof=1),np.mean((vals-mu)/sd)))
G=pd.DataFrame(rows,columns=['gene','fit','drug','p_t','sd','z'])
def bh(p):
 p=np.asarray(p); n=len(p); o=np.argsort(p); q=np.empty(n); q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]; return np.minimum(q,1)
G['q']=bh(G.p_t); G['call']=np.where(G.q<=.05,np.where(G.drug<0,'sensitiser','resistance'),'none'); G.to_csv('/task/output/rep_metrics.tsv',sep='\t',index=False)
X=R.merge(G,on='gene'); mods=X[X.expected_class!='none']; non=X[X.expected_class=='none']; print('t-test correct',sum(mods.call==mods.expected_class),'false none',sum(non.call!='none'),'calls',sum(G.call!='none'),'wrong',sum((mods.call!='none')&(mods.call!=mods.expected_class))); print(X.to_string(index=False))
