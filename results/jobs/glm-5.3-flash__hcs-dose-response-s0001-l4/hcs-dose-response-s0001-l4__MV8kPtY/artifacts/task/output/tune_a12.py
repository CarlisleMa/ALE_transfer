import numpy as np,pandas as pd,tifffile,multiprocessing as mp
from cellpose import models
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
M=None

def init():
 global M; M=models.Cellpose(gpu=False,model_type='nuclei')
def work(a):
 diam,cp,w=a
 x=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32); x=np.maximum((x-D)/f,0)
 masks,_,_,_=M.eval(x,diameter=diam,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=cp,min_size=15)
 regs=[r for r in regionprops(masks) if r.bbox[2]<320 and r.bbox[3]<320]
 areas=np.array([r.area for r in regs])
 pts=pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t')[['x_px','y_px']].values
 inside=[]
 for px,py in pts:
  xi=min(max(int(round(px)),0),319); yi=min(max(int(round(py)),0),319); inside.append(masks[yi,xi]>0)
 return diam,cp,w,len(regs),len(pts),sum(inside),np.percentile(areas,[0,10,50,90,100]).round(1).tolist(),int((areas>=30).sum())
if __name__=='__main__':
 jobs=[(d,cp,'A12') for cp in (-2,-1,0,0.5) for d in (6.5,7,7.5,8,9,10)]
 with mp.Pool(4,initializer=init) as pool:
  for r in pool.imap_unordered(work,jobs): print(r,flush=True)
