import numpy as np,pandas as pd
from scipy.optimize import least_squares
d=pd.read_csv('/task/output/explore_well_measurements.tsv',sep='\t')
d['row']=d.well.str[0]; d['col']=d.well.str[1:].astype(int)
qc=['C04','C05','E09','G06','G08']; u=d.query('well not in @qc').copy(); veh=u.query('role=="vehicle"')

def medpolish(df,val):
 tab=df.pivot(index='row',columns='col',values=val)
 grand=tab.stack().median(); re=pd.Series(0.,index=tab.index); ce=pd.Series(0.,index=tab.columns)
 for _ in range(200):
  re += (tab.sub(grand+ce,axis=1)).median(axis=1)-re
  grand=(tab.sub(re,axis=0).sub(ce)).stack().median()
  ce += (tab.sub(re,axis=0)).median(axis=0)-ce
  re += (tab.sub(grand+ce,axis=1)).median(axis=1)-re
  grand=(tab.sub(re,axis=0).sub(ce)).stack().median()
  ce += (tab.sub(re,axis=0)).median(axis=0)-ce
 return {f'{r}{c:02d}':grand+re[r]+ce[c] for r in tab.index for c in tab.columns}
def predmap(metric):
 return {
  'global':lambda g:np.repeat(veh[metric].mean(),len(g)),
  'column':lambda g:g.col.map(veh.groupby('col')[metric].mean()).values,
  'polish':lambda g:g.well.map(medpolish(veh,metric)).values,
 }
def fit(doses,responses):
 c=np.asarray(doses,float); y=np.asarray(responses,float); cmax=c.max()
 # bounds: log-ish EC50 through transform not needed with broad linear bounds
 def fun(p):
  b,logc,h=p
  ec=np.exp(logc); pred=b+(1-b)/(1+(c/ec)**h)
  return pred-y
 best=None
 for lc0 in np.linspace(np.log(c.min()/10),np.log(cmax*10),20):
  for b0 in (0,.05,.15,.3,.45):
   try:
    z=least_squares(fun,[b0,lc0,1],bounds=([0,np.log(c.min()/100),.05],[.999,np.log(cmax*1000),100]),max_nfev=10000)
    if best is None or z.cost<best.cost: best=z
   except Exception: pass
 if best is None:return None
 b,lc,h=best.x; ec=np.exp(lc)
 # solve response=.5: b+(1-b)/(1+(x/ec)^h)=.5 => x=ec*((1-b)/(1-2b?)) actually
 if b>=.5:return dict(bottom=b,ec50=ec,h=h,cross=np.inf,ok=False)
 cross=ec*((1-b)/(1-2*b))**(1/h)
 return dict(bottom=b,ec50=ec,h=h,cross=cross,ok=cross<=cmax)
metrics=['n_area40','r8_bgmode','r10_bgmode']
for metric in metrics:
 for norm,fn in predmap(metric).items():
  rows=[]
  for comp,g in u.query('role=="compound"').groupby('compound'):
   resp=g[metric]/fn(g); z=fit(g.dose_uM,resp)
   rows.append((comp,float(resp.min()),z['bottom'],z['h'],z['ec50'],z['cross'],'active' if z['ok'] else 'inactive'))
  print('\n###',metric,norm)
  print(pd.DataFrame(rows,columns=['compound','min_resp','bottom','h','ec50','cross','decision']).round(3).to_string(index=False))
