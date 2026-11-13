import numpy as np,pandas as pd,tifffile,multiprocessing as mp
from cellpose import models
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fl=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; fl.append(f/np.median(f))
M=None
def init():
 global M; M=models.CellposeModel(gpu=False,model_type='cyto3')
def work(a):
 diam,w=a
 ims=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  ims.append(np.maximum((raw-D)/fl[ch-1],0))
 masks,flows,styles=M.eval(np.stack(ims),channels=[2,1],diameter=diam,normalize=True,flow_threshold=.6,cellprob_threshold=0,min_size=30)
 regs=[r for r in regionprops(masks,intensity_image=ims[1]) if r.bbox[2]<320 and r.bbox[3]<320]
 # compare count to nuclei under final count and overlap
 nm,_,_=models.CellposeModel(gpu=False,model_type='nuclei').eval(ims[0],diameter=8.5,channels=[0,0],normalize=True,cellprob_threshold=-1,flow_threshold=2,min_size=15)
 nkeep=[r.label for r in regionprops(nm) if r.bbox[2]<320 and r.bbox[3]<320 and r.area>=40]
 vals=[]
 for lab in nkeep:
  nmask=nm==lab
  overlaps=np.bincount(masks[nmask].astype(int),minlength=int(masks.max())+1)
  overlaps[0]=0
  vals.append(int(overlaps.argmax()) if overlaps.max()>0 else 0)
 assigned=len(set(v for v in vals if v>0))
 areas=np.array([r.area for r in regs])
 return w,diam,'cyto_objects',len(regs),'nuclei_count',len(nkeep),'assigned_nuclei',assigned,'unique_cell_masks',len(set(vals)),'area_q',np.percentile(areas,[0,10,50,90,100]).round(1).tolist()
if __name__=='__main__':
 jobs=[(d,w) for d in (12,15,18) for w in ('A11','E05','A12','H12')]
 with mp.Pool(4,initializer=init) as pool:
  for r in pool.imap_unordered(work,jobs): print(*r,flush=True)
