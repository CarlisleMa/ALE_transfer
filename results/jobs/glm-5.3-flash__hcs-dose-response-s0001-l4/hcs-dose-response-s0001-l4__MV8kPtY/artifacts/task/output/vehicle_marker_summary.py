import numpy as np,pandas as pd,tifffile,glob,os
from scipy import ndimage as ndi
from skimage.measure import regionprops
from cellpose import models
from multiprocessing import Pool
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fl=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; fl.append(f/np.median(f))
M=None
def init(): global M; M=models.Cellpose(gpu=False,model_type='nuclei')
def work(w):
 ims=[]
 for ch in (1,2):
  raw=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32); ims.append(np.maximum((raw-D)/fl[ch-1],0))
 masks,_,_,_=M.eval(ims[0],diameter=10,channels=[0,0],normalize=True,cellprob_threshold=-1,flow_threshold=2,min_size=15)
 keep=[]
 for r in regionprops(masks):
  if r.bbox[2]<320 and r.bbox[3]<320 and r.area>=50: keep.append(r.label)
 nuc=np.isin(masks,keep); npx=keep
 if not keep: return dict(well=w,n=0,**{f'r{r}':np.nan for r in (8,10,12,14,16)},bg=np.nan)
 dist,inds=ndi.distance_transform_edt(~nuc,return_indices=True)
 nearest=np.where(nuc,masks,masks[inds[0],inds[1]])
 # robust extracellular background: histogram mode restricted to unassigned nonnuclear pixels
 bg=np.median(ims[1][~nuc])
 vals={}
 for rad in (8,10,12,14,16):
  cellvals=[]
  for lab in keep:
   m=(nearest==lab)&(~nuc)&(dist<=rad)
   cellvals.append(ims[1][m].mean()-bg)
  vals[f'r{rad}']=float(np.mean(cellvals))
 return dict(well=w,n=len(keep),bg=float(bg),**vals)
layout=pd.read_csv(base+'/plate_layout.tsv',sep='\t')
wells=layout.query('role=="vehicle" and well not in ["C04","C05","E09","G06","G08"]').well.tolist()
if __name__=='__main__':
 rows=[]
 with Pool(4,initializer=init) as pool:
  for r in pool.imap_unordered(work,wells): rows.append(r)
 df=pd.DataFrame(rows).merge(layout[['well','compound','dose_uM']],on='well')
 print(df.sort_values('well').to_string(index=False))
 print('\nmedians',df[[c for c in df if c.startswith('r')]].median().round(2).to_dict())
 df.to_csv('/task/output/vehicle_marker_radii.tsv',sep='\t',index=False)
