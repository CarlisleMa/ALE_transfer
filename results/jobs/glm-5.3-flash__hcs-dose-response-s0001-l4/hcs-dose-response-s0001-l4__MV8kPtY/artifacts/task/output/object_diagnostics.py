import glob,os,numpy as np,pandas as pd,tifffile
from scipy import ndimage as ndi
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f1=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f1/=np.median(f1)
f2=tifffile.imread(base+'/flatfield/flat_ch2.tif').astype(np.float32)-D; f2/=np.median(f2)
# fixed Cellpose settings, diameter 10
from cellpose import models
m=models.Cellpose(gpu=False,model_type='nuclei')
print('well point detected dist<=2 dist<=4 area q area Intensity Ch1 Ch2, centroid')
for p in sorted(glob.glob(base+'/annotations/*.tsv')):
 w=os.path.basename(p)[:-4]
 raw1=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32); x1=np.maximum((raw1-D)/f1,0)
 raw2=tifffile.imread(f'{base}/images/{w}_ch2.tif').astype(np.float32); x2=np.maximum((raw2-D)/f2,0)
 masks,_,_,_=m.eval(x1,diameter=10,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=0)
 pts=pd.read_csv(p,sep='\t')[['x_px','y_px']].values
 regs=regionprops(masks,intensity_image=x1)
 assigned={}
 stats=[]
 for r in regs:
  cy,cx=r.centroid
  ds=np.hypot(pts[:,0]-cx,pts[:,1]-cy)
  j=ds.argmin(); dist=ds[j]
  if dist<=4 and (r.label not in assigned or dist<assigned[r.label][1]):
   assigned[int(r.label)]=(int(j),float(dist))
  # local corrected marker: annulus 6..12 around nucleus, subtract edge/background median perhaps
  rr,cc=np.ogrid[:320,:320]
  ann=((rr-cy)**2+(cc-cx)**2>=36)&((rr-cy)**2+(cc-cx)**2<=144)
