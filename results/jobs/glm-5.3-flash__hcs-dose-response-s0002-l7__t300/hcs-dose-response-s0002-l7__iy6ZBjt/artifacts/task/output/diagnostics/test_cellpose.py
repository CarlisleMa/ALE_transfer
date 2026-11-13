import os, glob, time, numpy as np, pandas as pd, tifffile
from cellpose import models
from skimage import measure
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(np.float32)
F=tifffile.imread('/task/input/flatfield/flat_ch1.tif').astype(np.float32); F=np.maximum(F-D,0); F/=np.median(F)
wells=sorted(os.path.basename(p)[:-4] for p in glob.glob('/task/input/annotations/*.tsv'))
imgs=[]
for w in wells:
    imgs.append(np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(np.float32)-D,0)/F)
t=time.time(); model=models.CellposeModel(gpu=False,model_type='nuclei')
masks,flows,styles=model.eval(imgs,channels=[[0,0]]*len(imgs),diameter=14.0,flow_threshold=0.6,cellprob_threshold=0.0,batch_size=8)
print('elapsed',time.time()-t,flush=True)
np.savez_compressed('/task/output/diagnostics/cellpose_d14_masks.npz',**{w:m for w,m in zip(wells,masks)})
for w,mask in zip(wells,masks):
 manual=pd.read_csv(f'/task/input/annotations/{w}.tsv',sep='\t')
 regs=measure.regionprops(mask)
 pts=manual[['x_px','y_px']].to_numpy()
 # assign each annotated point to containing mask if possible; otherwise nearest mask centroid
 cents=np.array([r.centroid for r in regs])
 assigned=[]; unmatched_pts=[]
 for x,y in pts:
    lab=mask[int(round(y)),int(round(x))]
    if lab: assigned.append(lab-1)
    elif len(regs):
       j=np.argmin((cents[:,0]-y)**2+(cents[:,1]-x)**2)
       assigned.append(j); unmatched_pts.append((x,y))
 used=set(assigned)
 fp=[]
 for j,r in enumerate(regs):
    if j not in used: fp.append((int(r.area),round(float(r.minor_axis_length),1),round(float(r.major_axis_length),1),round(float(r.centroid[1]),1),round(float(r.centroid[0]),1)))
 print(w,'manual',len(manual),'raw',len(regs),'unmatched_pts',len(unmatched_pts),'unmatched_pts_dist',[(p,round(float(np.min((cents[:,0]-p[1])**2+(cents[:,1]-p[0])**2)),1)) for p in unmatched_pts],flush=True)
 print(' FP area/minor/major/x/y',sorted(fp),flush=True)
