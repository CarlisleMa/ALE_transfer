import numpy as np,pandas as pd,tifffile,multiprocessing as mp
from cellpose import models
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
M=None
def init():
 global M; M=models.Cellpose(gpu=False,model_type='nuclei')
def work(a):
 diam,cp,ft,w=a
 x=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32); x=np.maximum((x-D)/f,0)
 masks,_,_,_=M.eval(x,diameter=diam,channels=[0,0],normalize=True,cellprob_threshold=cp,flow_threshold=ft,min_size=15)
 regs=[r for r in regionprops(masks) if r.bbox[2]<320 and r.bbox[3]<320]
 areas=np.array([r.area for r in regs])
 pts=pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t')[['x_px','y_px']].values
 inside=[]
 for px,py in pts:
  xi=min(max(int(round(px)),0),319); yi=min(max(int(round(py)),0),319); inside.append(bool(masks[yi,xi]>0))
 vals=[len(regs)]+[int((areas>=z).sum()) for z in (25,30,40,50,55,60,70,80)]
 return diam,cp,ft,w,'manual',len(pts),'raw_and_filters',vals,'points_inside',sum(inside)
if __name__=='__main__':
 wells=['A06','A11','A12','C12','D10','E05','G07','H12']
 jobs=[(d,cp,ft,w) for d in (8.5,10) for cp in (-1,0) for ft in (.6,2) for w in wells]
 with mp.Pool(4,initializer=init) as pool:
  for r in pool.imap_unordered(work,jobs): print(*r,flush=True)
