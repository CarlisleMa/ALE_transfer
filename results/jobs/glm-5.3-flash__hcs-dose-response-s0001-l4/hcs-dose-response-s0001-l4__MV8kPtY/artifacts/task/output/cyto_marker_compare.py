import numpy as np,pandas as pd,tifffile,multiprocessing as mp
from scipy import ndimage as ndi
from skimage.measure import regionprops
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
FL=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; FL.append(f/np.median(f))
NC=None; CC=None
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
 keep=[r.label for r in regionprops(nm) if r.bbox[2]<320 and r.bbox[3]<320 and r.area>=40]
 nuc=np.isin(nm,keep)
 cm,_,_=CC.eval(np.stack(ims),channels=[2,1],diameter=15,normalize=True,flow_threshold=.6,cellprob_threshold=0,min_size=30)
 # assign each kept nucleus to maximum-overlap cyto3 mask
 vals=[]; used_cells=set(); unassigned=0
 for lab in keep:
  nmask=nm==lab
  counts=np.bincount(cm[nmask].astype(int),minlength=int(cm.max())+1); counts[0]=0
  target=int(counts.argmax())
  if counts[target]==0: unassigned+=1; continue
  used_cells.add(target)
  cell=(cm==target)&(~nuc)
  vals.append((w,lab,counts[target],cell.sum(),float(ims[1][cell].mean()) if cell.any() else np.nan))
 # territories
 dist,inds=ndi.distance_transform_edt(~nuc,return_indices=True)
 nearest=np.where(nuc,nm,nm[inds[0],inds[1]])
 bg=hist_mode(ims[1])
 out={'well':w,'n_counted':len(keep),'cyto_cells':int(cm.max()),'assigned_nuclei':len(keep)-unassigned,'unique_cyto_used':len(used_cells),'bg':bg}
 cy=np.array([v[4] for v in vals]) if vals else np.array([])
 # Build territory means
 terr={}
 for rad in (8,10,12):
  per=[]
  for lab in keep:
   m=(nearest==lab)&(~nuc)&(dist<=rad)
   if m.any(): per.append(float(ims[1][m].mean()))
  terr[f'r{rad}']=np.array(per)
 out['cyto_mean_bg']=float(np.mean(cy)-bg) if len(cy) else np.nan
 for rad,a in terr.items(): out[f'{rad}_mean_bg']=float(np.mean(a)-bg) if len(a) else np.nan
 # per-cell correlation between cyto and territories on shared nuclei
 shared=[]
 for i,v in enumerate(vals):
  lab=v[1]
  # recompute territory vector
  row=[v[4]]
  for rad,a in terr.items(): row.append(a[i])
  shared.append(row)
 arr=np.asarray(shared)
 if len(arr)>2 and np.std(arr[:,0])>0:
  out['corr_cyto_r8']=float(np.corrcoef(arr[:,0],arr[:,1])[0,1]); out['corr_cyto_r10']=float(np.corrcoef(arr[:,0],arr[:,2])[0,1])
 return out
if __name__=='__main__':
 wells=['A03','E01','A05','H03','C12','G12','F12','E12']
 with mp.Pool(4,initializer=init) as pool:
  for x in pool.imap_unordered(work,wells): print(x,flush=True)
