import numpy as np,pandas as pd,tifffile
from cellpose import models
from scipy import ndimage as ndi
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32); f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
x=tifffile.imread(base+'/images/A12_ch1.tif').astype(np.float32); x=np.maximum((x-D)/f,0)
m=models.CellposeModel(gpu=False,model_type='nuclei')
_,flows,_=m.eval(x,diameter=10,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=-2)
p=flows[2]
pts=pd.read_csv(base+'/annotations/A12.tsv',sep='\t')[['x_px','y_px']].values
print('cellprob global q',np.percentile(p,[0,1,5,10,50,90,95,99,100]))
for i,(px,py) in enumerate(pts):
 xi,yi=int(round(px)),int(round(py))
 print(i,px,py,'prob3x3',np.round(p[yi-1:yi+2,xi-1:xi+2],3),'patch max',round(p[max(0,yi-8):yi+9,max(0,xi-8):xi+9].max(),3),'ch1 patch q/max',np.percentile(x[max(0,yi-8):yi+9,max(0,xi-8):xi+9],[1,50,90,99,100]).round(1))
