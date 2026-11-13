import numpy as np,pandas as pd,tifffile,multiprocessing as mp
from scipy import ndimage as ndi
from skimage.measure import regionprops
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
FL=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; FL.append(f/np.median(f))
NC=CC=None
def init():
 global NC,CC; NC=models.CellposeModel(gpu=False,model_type='nuclei'); CC=models.CellposeModel(gpu=False,model_type='cyto3')
def hist_mode(x):
 lo=max(0,np.floor(np.percentile(x,.5))); hi=np.percentile(x,99.5)
 bins=np.arange(lo,max(lo+3,hi+2),2.0)
 h,e=np.histogram(x,bins=bins); h=ndi.uniform_filter1d(h.astype(float),3)
 return float((e[np.argmax(h)]+e[np.argmax(h)+1])/2)
def work(w):
 ims=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  ims.append(np.maximum((raw-D)/FL[ch-1],0))
 nm,_,_=NC.eval(ims[0],diameter=8.5,channels=[0,0],normalize=True,cellprob_threshold=-1,flow_threshold=2,min_size=15)
 keep=[r.label for r in regionprops(nm) if r.bbox[2]<320 and r.bbox[3]<320 and r.area>=50]
 nuc=np.isin(nm,keep)
 # nearest-nucleus labels for all pixels
 dist,inds=ndi.distance_transform_edt(~nuc,return_indices=True)
 nearest=np.where(nuc,nm,nm[inds[0],inds[1]])
 # cyto3 foreground and partition by nearest counted nucleus
 cm,_,_=CC.eval(np.stack(ims),channels=[2,1],diameter=15,normalize=True,flow_threshold=.6,cellprob_threshold=0,min_size=30)
 fg=cm>0
 cyto_partition=fg&(~nuc)
 # Remove pixels whose nearest nucleus is not counted (possible in rare gaps)
 cyto_partition&=np.isin(nearest,keep)
 bg=hist_mode(ims[1])
 vals={k:[] for k in ['r8','r10','cytopart']}
 no_cyto=[]
 for lab in keep:
  nucm=nm==lab
  m8=(nearest==lab)&(~nuc)&(dist<=8)
  m10=(nearest==lab)&(~nuc)&(dist<=10)
  mp=(nearest==lab)&cyto_partition
  if not mp.any():
   # fallback to marker/cyto foreground within 12 pixels
   mp=(nearest==lab)&(~nuc)&(dist<=12)&fg
  if not mp.any(): no_cyto.append(int(lab)); mp=m8
  vals['r8'].append(float(ims[1][m8].mean()) if m8.any() else np.nan)
  vals['r10'].append(float(ims[1][m10].mean()) if m10.any() else np.nan)
  vals['cytopart'].append(float(ims[1][mp].mean()))
 return dict(well=w,n=len(keep),bg=bg,no_cyto=len(no_cyto),**{f'{k}_mean_bg_sub':float(np.nanmean(v)-bg) for k,v in vals.items()})
if __name__=='__main__':
 wells=['A03','E01','A05','H03','C12','G12','F12','E12','A11','E05']
 with mp.Pool(4,initializer=init) as pool:
  for x in pool.imap_unordered(work,wells): print(x,flush=True)
