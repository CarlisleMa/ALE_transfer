import glob,os,numpy as np,pandas as pd,tifffile,time,multiprocessing as mp
from scipy.spatial import ConvexHull
from skimage.measure import regionprops
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fl=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; fl.append(f/np.median(f))

def min_feret(mask):
 ys,xs=np.nonzero(mask); pts=np.c_[xs,ys].astype(float)
 h=pts[ConvexHull(pts).vertices]
 vals=[]
 for i in range(len(h)):
  p=h[i]; q=h[(i+1)%len(h)]; v=q-p; L=np.hypot(*v)
  if L>0: vals.append(np.max((h-p)@(v/L)))
 return min(vals) if vals else 0

def work(args):
 diam,w,pth=args
 t=time.time()
 x=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  x.append(np.maximum((raw-D)/fl[ch-1],0))
 model=models.Cellpose(gpu=False,model_type='nuclei')
 masks,_,_,_=model.eval(x[0],diameter=diam,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=0)
 pts=pd.read_csv(pth,sep='\t')[['x_px','y_px']].values
 regs=regionprops(masks,intensity_image=x[0])
 # one-to-one nearest match; candidates by distance <=8
 pairs=[]
 for ri,r in enumerate(regs):
  cy,cx=r.centroid
  for pi,(px,py) in enumerate(pts):
   d=np.hypot(px-cx,py-cy)
   if d<=8: pairs.append((d,ri,pi))
 usedo=set(); usedp=set()
 for d,ri,pi in sorted(pairs):
  if ri not in usedo and pi not in usedp: usedo.add(ri); usedp.add(pi)
 out=[]
 for ri,r in enumerate(regs):
  cy,cx=r.centroid; sl=masks[r.slice]==r.label
  rr,cc=np.ogrid[:320,:320]; ring=((rr-cy)**2+(cc-cx)**2>=30)&((rr-cy)**2+(cc-cx)**2<=160)
  j=next((pi for d,oo,pi in sorted(pairs) if oo==ri and pi not in usedp),None) # not meaningful for used
  dist=np.min(np.hypot(pts[:,0]-cx,pts[:,1]-cy))
  matched=ri in usedo
  out.append({'diameter':diam,'well':w,'label':int(r.label),'matched':matched,'area':int(r.area),'equiv_diam':float(2*np.sqrt(r.area/np.pi)),
   'bbox_min':float(min(r.bbox[2]-r.bbox[0],r.bbox[3]-r.bbox[1])),'axis_minor':float(getattr(r,'axis_minor_length',np.nan)),
   'min_feret':float(min_feret(sl)),'dist':float(dist),'ch1_mean':float(x[0][r.slice][sl].mean()),
   'ch2_nucleus':float(x[1][r.slice][sl].mean()),'ch2_ring5_12':float(x[1][ring].mean()),'cy':float(cy),'cx':float(cx)})
 return out

def run(diam,w,pth):
 try: return work((diam,w,pth))
 except Exception as e: return [{'diameter':diam,'well':w,'error':str(e)}]
if __name__=='__main__':
 jobs=[]
 for diam in (8.5,10):
  for pth in sorted(glob.glob(base+'/annotations/*.tsv')):
   w=os.path.basename(pth)[:-4]; jobs.append((diam,w,pth))
 allrows=[]
 with mp.Pool(4) as pool:
  for rows in pool.imap_unordered(work,jobs):
   allrows.extend(rows)
   if rows and 'error' in rows[0]: print('ERROR',rows,flush=True)
 df=pd.DataFrame(allrows)
 df.to_csv('/task/output/nuclear_objects.tsv',sep='\t',index=False)
 for (d,w),g in df.groupby(['diameter','well']):
  manual=len(pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t'))
  print(d,w,'manual',manual,'objects',len(g),'matched',g.matched.sum(),'unmatched_obj',(~g.matched).sum(),'missed_points',manual-g.matched.sum(),flush=True)
  for matched in (True,False):
   a=g[g.matched==matched]
   if len(a): print(' ', 'matched' if matched else 'extra',len(a),'area q10/50/90',a.area.quantile([.1,.5,.9]).round(1).tolist(),'eqdiam',a.equiv_diam.quantile([.1,.5,.9]).round(1).tolist(),'minor',a.axis_minor.quantile([.1,.5,.9]).round(1).tolist(),'feret',a.min_feret.quantile([.1,.5,.9]).round(1).tolist(),'ch2ring',a.ch2_ring5_12.quantile([.1,.5,.9]).round(1).tolist())
