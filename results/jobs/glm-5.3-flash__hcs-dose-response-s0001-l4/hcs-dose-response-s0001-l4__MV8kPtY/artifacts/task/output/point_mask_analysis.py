import glob,os,numpy as np,pandas as pd,tifffile
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
for p in sorted(glob.glob(base+'/annotations/*.tsv')):
 w=os.path.basename(p)[:-4]; z=np.load(f'/task/output/masks_d10_{w}.npz'); m=z['masks']; x=z['corrected']
 pts=pd.read_csv(p,sep='\t')[['x_px','y_px']].values
 regs={r.label:r for r in regionprops(m,intensity_image=x)}
 labels=[]
 details=[]
 for px,py in pts:
  xi,yi=int(round(px)),int(round(py)); xi=min(max(xi,0),319); yi=min(max(yi,0),319)
  lab=int(m[yi,xi]); labels.append(lab)
  if lab==0:
   # nearest centroid
   ds=[(np.hypot(px-r.centroid[1],py-r.centroid[0]),r.label) for r in regs.values()]
   nd,lab=min(ds) if ds else (np.nan,0)
  details.append((lab,px,py))
 unique,cnt=np.unique([a for a in labels if a>0],return_counts=True)
 print('\n',w,'manual',len(pts),'point-covered masks',len(unique),'multiplicities',dict(zip(unique.astype(int),cnt.astype(int))))
 for lab,n in zip(unique,cnt):
  r=regs[lab]; print(' mask',int(lab),'npts',int(n),'area',r.area,'bbox',r.bbox,'minor',round(getattr(r,'axis_minor_length',np.nan),1),'mean',round(x[r.slice][m[r.slice]==lab].mean(),1))
 # unmatched points
 for (lab,px,py) in details:
  if lab==0: print(' missed point',px,py)
