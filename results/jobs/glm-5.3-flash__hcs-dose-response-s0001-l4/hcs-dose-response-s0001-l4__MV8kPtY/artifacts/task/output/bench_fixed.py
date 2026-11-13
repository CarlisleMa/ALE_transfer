import glob, os, numpy as np, pandas as pd, tifffile, time, multiprocessing as mp
from cellpose import models
from skimage.measure import regionprops
base='/task/input'
DARK=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
FLAT=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-DARK
FLAT/=np.median(FLAT)
def corrected(w):
 x=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32)
 return np.maximum((x-DARK)/FLAT,0)
def frame(r):
 return r.bbox[2]<320 and r.bbox[3]<320

def work(args):
 diam,w,pth=args
 m=models.Cellpose(gpu=False,model_type='nuclei')
 im=corrected(w)
 t=time.time(); masks,_,_,dd=m.eval(im,diameter=diam,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=0)
 regs=[r for r in regionprops(masks) if frame(r)]
 manual=len(pd.read_csv(pth,sep='\t'))
 stats=[]
 for r in regs:
  # robust min Feret: minimum distance between opposite supporting convex hull points
  from scipy.spatial import ConvexHull
  ys,xs=np.nonzero(masks[r.slice]==r.label); pts=np.c_[xs,ys]
  h=pts[ConvexHull(pts).vertices]
  d=[]
  for i in range(len(h)):
   p=h[i]; q=h[(i+1)%len(h)]
   v=q-p; L=np.hypot(*v)
   if L:
    d.append(np.max(np.abs((h-p)@(v/L))))
  stats.append(min(d) if d else 0)
 return w,manual,len(regs),np.percentile(stats,[0,5,10,25,50,75,95,100]).round(1).tolist(),int((time.time()-t)*1000)

if __name__=='__main__':
 jobs=[]
 for diam in (9.2,10,11):
  for pth in sorted(glob.glob(base+'/annotations/*.tsv')):
   w=os.path.basename(pth)[:-4]; jobs.append((diam,w,pth))
 with mp.Pool(8) as pool:
  for x in pool.imap_unordered(work,jobs):
   print(x,flush=True)
