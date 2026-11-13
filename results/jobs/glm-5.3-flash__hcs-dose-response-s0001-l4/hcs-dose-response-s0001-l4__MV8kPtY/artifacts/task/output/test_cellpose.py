import glob, os, numpy as np, pandas as pd, tifffile, time
from cellpose import models
from skimage.measure import label, regionprops
from scipy.ndimage import gaussian_laplace
base='/task/input'
DARK=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
flats={ch:(tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-DARK) for ch in (1,2)}
def corr(im,ch):
    x=(im.astype(np.float32)-DARK)/flats[ch]
    return x*np.median(flats[ch])
def longest_shortest(mask):
    # Approximate minimum Feret via area and major-axis-derived mean width; use robust min-distance boundary method later if needed.
    from skimage.measure import regionprops
    return None
def passes_frame(r):
    minr,minc,maxr,maxc=r.bbox
    return not (maxr==320 or maxc==320) # bottom/right exclusion; top/left inclusion
model=models.Cellpose(gpu=False,model_type='nuclei')
for diam in [0,10,12,15,18,22]:
 print('\nDIAM',diam,flush=True)
 for p in sorted(glob.glob(base+'/annotations/*.tsv')):
  w=os.path.basename(p)[:-4]
  im=corr(tifffile.imread(f'{base}/images/{w}_ch1.tif'),1)
  t=time.time(); masks,flows,styles,diams=model.eval(im,diameter=diam,channels=[0,0],normalize=True,flow_threshold=0.6,cellprob_threshold=0.0)
  regs=[r for r in regionprops(masks) if passes_frame(r)]
  # area and bbox width filters (rough first pass)
  keep=[r for r in regs if r.area>=20]
  manual=len(pd.read_csv(p,sep='\t'))
  print(w,'manual',manual,'raw',len(regs),'area20',len(keep),'msec',int((time.time()-t)*1000),'diam_used',diams,flush=True)
