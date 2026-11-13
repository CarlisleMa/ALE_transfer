import os,glob,numpy as np,pandas as pd,tifffile,multiprocessing as mp
from scipy import ndimage as ndi, spatial
from skimage.measure import regionprops
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
FL=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; FL.append(f/np.median(f))
M=None
def init():
 global M; M=models.Cellpose(gpu=False,model_type='nuclei')
def min_width(mask):
 ys,xs=np.nonzero(mask); pts=np.c_[xs,ys].astype(float)
 if len(pts)<3: return float(min(np.ptp(xs),np.ptp(ys))+1)
 h=pts[spatial.ConvexHull(pts).vertices]; vals=[]
 for i in range(len(h)):
  p=h[i]; q=h[(i+1)%len(h)]; v=q-p; L=np.hypot(*v)
  if L:
   n=np.array([-v[1],v[0]])/L; vals.append(np.ptp((h-p)@n))
 return min(vals) if vals else 0.0
def hist_mode(x):
 lo=np.floor(np.percentile(x,0.5)); hi=np.percentile(x,99.5)
 bins=np.arange(max(0,lo),max(lo+3,hi+2),2.0)
 h,e=np.histogram(x,bins=bins); h=ndi.uniform_filter1d(h.astype(float),3)
 k=np.argmax(h)
 return float((e[k]+e[k+1])/2)
def work(w):
 ims=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  ims.append(np.maximum((raw-D)/FL[ch-1],0))
 masks,_,_,_=M.eval(ims[0],diameter=8.5,channels=[0,0],normalize=True,cellprob_threshold=-1,flow_threshold=2,min_size=15)
 objs=[]
 for r in regionprops(masks,intensity_image=ims[0]):
  if r.bbox[2]==320 or r.bbox[3]==320: continue
  sl=masks[r.slice]==r.label
  objs.append((int(r.label),int(r.area),float(min_width(sl)),float(getattr(r,'axis_minor_length',np.nan)),float(ims[0][r.slice][sl].mean()),float(r.bbox[0]),float(r.bbox[1]),float(r.bbox[2]),float(r.bbox[3])))
 o=pd.DataFrame(objs,columns=['label','area','min_width','axis_minor','ch1','r0','c0','r1','c1'])
 keep=o[(o.area>=40)|(o.min_width>=11/1.3)].label.astype(int).to_numpy()
 nucmask=np.isin(masks,keep)
 bgmode=hist_mode(ims[1]); bgp10=float(np.percentile(ims[1],10)); bgnon=float(np.median(ims[1][~nucmask])) if (~nucmask).any() else bgmode
 vals={'well':w,'focus_ten':float(np.mean(sum(g*g for g in np.gradient(ndi.gaussian_filter(ims[0],0.5))))),'n_frame':len(o),'n_area40':int((o.area>=40).sum()),'n_area50':int((o.area>=50).sum()),'n_area_or_width':len(keep),'n_width846':int((o.min_width>=11/1.3).sum()),'bg_mode':bgmode,'bg_p10':bgp10,'bg_nonnuc':bgnon}
 if len(keep):
  dist,inds=ndi.distance_transform_edt(~nucmask,return_indices=True)
  nearest=np.where(nucmask,masks,masks[inds[0],inds[1]])
  for radius in (8,10,12):
   per=[]
   for lab in keep:
    m=(nearest==lab)&(~nucmask)&(dist<=radius)
    if m.any(): per.append(float(ims[1][m].mean()))
   vals[f'r{radius}_raw']=float(np.mean(per)); vals[f'r{radius}_bgmode']=float(np.mean(per)-bgmode); vals[f'r{radius}_bgp10']=float(np.mean(per)-bgp10)
 else:
  for radius in (8,10,12):
   vals[f'r{radius}_raw']=np.nan; vals[f'r{radius}_bgmode']=np.nan; vals[f'r{radius}_bgp10']=np.nan
 return vals
if __name__=='__main__':
 layout=pd.read_csv(base+'/plate_layout.tsv',sep='\t'); wells=layout.well.tolist()
 rows=[]
 with mp.Pool(4,initializer=init) as pool:
  for i,r in enumerate(pool.imap_unordered(work,wells),1):
   rows.append(r)
   if i%12==0 or i==len(wells): print('done',i,flush=True)
 df=pd.DataFrame(rows).merge(layout,on='well').sort_values('well')
 df.to_csv('/task/output/explore_well_measurements.tsv',sep='\t',index=False)
 print(df[['well','compound','dose_uM','role','n_area40','n_area50','n_area_or_width','r8_bgmode','r10_bgmode','r12_bgmode']].to_string(index=False))
