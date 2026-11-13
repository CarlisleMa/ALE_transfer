import glob,os,numpy as np,pandas as pd,tifffile,multiprocessing as mp
from scipy.spatial import ConvexHull
from skimage.measure import regionprops
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fl=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; fl.append(f/np.median(f))
M=None
def init():
 global M; M=models.Cellpose(gpu=False,model_type='nuclei')
def min_width(mask):
 ys,xs=np.nonzero(mask); pts=np.c_[xs,ys].astype(float)
 h=pts[ConvexHull(pts).vertices]
 vals=[]
 for i in range(len(h)):
  p=h[i]; q=h[(i+1)%len(h)]; v=q-p; L=np.hypot(*v)
  if L>0:
   n=np.array([-v[1],v[0]])/L
   vals.append(np.ptp((h-p)@n))
 return min(vals) if vals else 0
def work(w):
 x=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  x.append(np.maximum((raw-D)/fl[ch-1],0))
 masks,_,_,_=M.eval(x[0],diameter=10,channels=[0,0],normalize=True,cellprob_threshold=-1,flow_threshold=2,min_size=15)
 np.savez_compressed(f'/task/output/final_masks_{w}.npz',masks=masks.astype(np.int32),ch1=x[0],ch2=x[1])
 pts=pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t')[['x_px','y_px']].values
 regs=[r for r in regionprops(masks,intensity_image=x[0]) if r.bbox[2]<320 and r.bbox[3]<320]
 pairs=[]
 for ri,r in enumerate(regs):
  cy,cx=r.centroid
  for pi,(px,py) in enumerate(pts):
   d=np.hypot(px-cx,py-cy)
   if d<=9: pairs.append((d,ri,pi))
 usedo=set(); usedp=set()
 for d,ri,pi in sorted(pairs):
  if ri not in usedo and pi not in usedp: usedo.add(ri); usedp.add(pi)
 rows=[]
 for ri,r in enumerate(regs):
  sl=masks[r.slice]==r.label; cy,cx=r.centroid
  rr,cc=np.ogrid[:320,:320]; dist2=(rr-cy)**2+(cc-cx)**2
  rings=[]
  for lo,hi in ((0,3),(3,6),(6,9),(9,12),(12,16),(16,22),(22,30)):
   rings.append(float(x[1][(dist2>=lo*lo)&(dist2<hi*hi)].mean()))
  rows.append({'well':w,'label':int(r.label),'matched':ri in usedo,'area':int(r.area),'min_width':float(min_width(sl)),'equiv_diam':float(2*np.sqrt(r.area/np.pi)),'axis_minor':float(getattr(r,'axis_minor_length',np.nan)),'ch1':float(x[0][r.slice][sl].mean()),'ch2_nucleus':float(x[1][r.slice][sl].mean()),**{f'ring_{lo}_{hi}':v for (lo,hi),v in zip(((0,3),(3,6),(6,9),(9,12),(12,16),(16,22),(22,30)),rings)}})
 return rows
if __name__=='__main__':
 wells=[os.path.basename(p)[:-4] for p in sorted(glob.glob(base+'/annotations/*.tsv'))]
 rows=[]
 with mp.Pool(4,initializer=init) as pool:
  for rr in pool.imap_unordered(work,wells): rows.extend(rr)
 df=pd.DataFrame(rows); df.to_csv('/task/output/final_object_diagnostics.tsv',sep='\t',index=False)
 for w,g in df.groupby('well'):
  manual=len(pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t'))
  print('\n',w,'manual',manual,'frame objects',len(g),'matched',g.matched.sum(),'extra',(~g.matched).sum())
  for status,h in g.groupby('matched'):
   print(status,'n',len(h),'area q',h.area.quantile([0,.05,.1,.25,.5,.75,.9,.95,1]).round(1).tolist(),'min_width q',h.min_width.quantile([0,.05,.1,.25,.5,.75,.9,.95,1]).round(1).tolist())
