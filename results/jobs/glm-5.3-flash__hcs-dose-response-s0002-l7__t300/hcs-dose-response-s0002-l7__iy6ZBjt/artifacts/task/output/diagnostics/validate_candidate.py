import glob,os,numpy as np,pandas as pd,tifffile
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(np.float32)
F=tifffile.imread('/task/input/flatfield/flat_ch1.tif').astype(np.float32); F=np.maximum(F-D,0); F/=np.median(F)

def segment(im):
 sm=ndi.gaussian_filter(im,0.75)
 pts=feature.peak_local_max(sm,min_distance=3,threshold_abs=900,exclude_border=False)
 pts=pts[~((pts[:,0]>=319)|(pts[:,1]>=319))]
 fg=ndi.binary_fill_holes(morphology.remove_small_objects(sm>filters.threshold_otsu(sm),10))
 markers=np.zeros(im.shape,np.int32)
 for i,(r,c) in enumerate(pts,1): markers[r,c]=i
 lab=segmentation.watershed(-sm,markers,mask=fg)
 keep=[]; stats=[]
 for r in measure.regionprops(lab):
  minr,minc,maxr,maxc=r.bbox
  if maxr>=320 or maxc>=320: continue
  if r.area>=35:
   keep.append(r.label)
   stats.append((r.label,r.area,r.minor_axis_length,r.major_axis_length,minr,minc,maxr,maxc))
 return lab,np.isin(lab,keep),pd.DataFrame(stats,columns=['label','area','minor','major','r0','c0','r1','c1'])

rows=[]
for p in sorted(glob.glob('/task/input/annotations/*.tsv')):
 w=os.path.basename(p)[:-4]
 im=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(np.float32)-D,0)/F
 lab,keep,stats=segment(im)
 manual=len(pd.read_csv(p,sep='\t')); count=len(stats); tol=max(3,.1*manual)
 rows.append(dict(well=w,manual=manual,pipeline=count,delta=count-manual,tol=tol,pass_validation=abs(count-manual)<=tol,
                  n_raw=int(lab.max()),area_p10=np.percentile(stats.area,10),area_p50=np.percentile(stats.area,50),area_p90=np.percentile(stats.area,90)))
d=pd.DataFrame(rows)
print(d.to_string(index=False,float_format=lambda v:f'{v:.3g}'))
print('\nALL VALIDATED:',bool(d.pass_validation.all()))
