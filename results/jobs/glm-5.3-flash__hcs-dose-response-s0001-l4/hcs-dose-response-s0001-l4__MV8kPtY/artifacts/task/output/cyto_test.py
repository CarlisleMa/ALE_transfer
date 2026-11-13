import numpy as np, pandas as pd, tifffile
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
fl=[]
for ch in (1,2):
 f=tifffile.imread(f'{base}/flatfield/flat_ch{ch}.tif').astype(np.float32)-D; f/=np.median(f); fl.append(f)
for w in ['A11','E05','A12','D10','H12']:
 ims=[]
 for ch in (1,2):
  x=tifffile.imread(f'{base}/images/{w}_ch{ch}.tif').astype(np.float32)
  ims.append(np.maximum((x-D)/fl[ch-1],0))
 img=np.stack(ims) # channels first
 m=models.CellposeModel(gpu=False,model_type='cyto3')
 masks,flows,styles=m.eval(img,channels=[1,2],diameter=18,normalize=True,flow_threshold=.6,cellprob_threshold=0)
 print(w,'cyto objects',masks.max(),'nuclear manual',len(pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t')) if w in pd.read_csv('/task/input/plate_layout.tsv',sep='\t').query('role=="compound"').well.values else '?','mask areas q',np.percentile(np.bincount(masks.ravel())[1:],[0,10,25,50,75,90,100]).round(1))
