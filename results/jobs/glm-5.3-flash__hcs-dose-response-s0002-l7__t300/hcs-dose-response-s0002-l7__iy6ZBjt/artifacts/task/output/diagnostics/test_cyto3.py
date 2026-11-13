import time, numpy as np, pandas as pd, tifffile
from cellpose import models
from skimage import measure
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
Fl={c:tifffile.imread(f'/task/input/flatfield/flat_ch{c}.tif').astype(float) for c in (1,2)}
for c in Fl: Fl[c]=np.maximum(Fl[c]-D,0)/np.median(Fl[c]-D)
model=models.CellposeModel(gpu=False,model_type='cyto3')
for w in ['A03','E02','B12']:
 im2=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch2.tif').astype(float)-D,0)/Fl[2]
 im1=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/Fl[1]
 t=time.time()
 masks,flows,styles=model.eval(np.stack([im2,im1]),channels=[1,2],channel_axis=0,diameter=20,flow_threshold=.6,cellprob_threshold=0)
 regs=measure.regionprops(masks)
 bg=np.median(im2[masks==0]) if np.any(masks==0) else np.nan
 means=np.array([float(im2[masks==r.label].mean()) for r in regs]) if regs else np.array([])
 print(w,'seconds',round(time.time()-t,1),'n_cells',len(regs),'areas',np.percentile([r.area for r in regs],[0,10,50,90,100]).round(0),'bg',round(bg,1),'cell mean median',round(float(np.median(means)),1) if len(means) else None,'corrected median',round(float(np.median(means)-bg),1) if len(means) else None,flush=True)
 np.save(f'/task/output/diagnostics/cyto3_{w}.npy',masks)
