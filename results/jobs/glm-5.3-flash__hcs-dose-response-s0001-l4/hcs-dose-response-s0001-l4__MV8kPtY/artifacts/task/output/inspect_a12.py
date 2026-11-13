import numpy as np,pandas as pd,tifffile
from scipy import ndimage as ndi
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fs=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; fs.append(f/np.median(f))
x=[]
for ch in (1,2):
 r=tifffile.imread(f'{base}/images/A12_ch{ch}.tif').astype(np.float32); x.append(np.maximum((r-D)/fs[ch-1],0))
pts=pd.read_csv(base+'/annotations/A12.tsv',sep='\t')[['x_px','y_px']].values
chars=' .:-=+*#%@'
for i,(px,py) in enumerate(pts):
 print('\nPOINT',i,'x y',px,py)
 for ch in (0,1):
  xi,yi=int(round(px)),int(round(py)); lo=max(0,xi-15);hi=min(320,xi+16); loy=max(0,yi-15);hiy=min(320,yi+16)
  a=x[ch][loy:hiy,lo:hi]
  print(' ch',ch+1,'3x3',np.round(x[ch][yi-1:yi+2,xi-1:xi+2],1),'patch q/min/max/mean/std',np.percentile(a,[1,10,50,90,99]).round(1),round(float(a.min()),1),round(float(a.max()),1),round(float(a.mean()),1),round(float(a.std()),1))
  # downsample patch to ~15x15, normalize 1-99 percentile
  b=a[::2,::2]; z=(b-np.percentile(b,2))/(np.percentile(b,98)-np.percentile(b,2)+1e-6)
  for row in z: print(''.join(chars[min(9,max(0,int(v*9)))] for v in row))
