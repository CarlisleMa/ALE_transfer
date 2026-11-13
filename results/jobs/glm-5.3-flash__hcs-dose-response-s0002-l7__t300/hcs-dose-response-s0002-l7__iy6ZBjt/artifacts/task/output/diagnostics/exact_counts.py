import numpy as np,pandas as pd,tifffile,glob,os
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
F=tifffile.imread('/task/input/flatfield/flat_ch1.tif').astype(float); F=np.maximum(F-D,0); F/=np.median(F)
layout=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')

def count(im):
 sm=ndi.gaussian_filter(im,.75)
 pts=feature.peak_local_max(sm,min_distance=3,threshold_abs=900,exclude_border=False)
 pts=pts[~((pts[:,0]>=319)|(pts[:,1]>=319))]
 fg=ndi.binary_fill_holes(morphology.remove_small_objects(sm>filters.threshold_otsu(sm),10))
 markers=np.zeros(im.shape,np.int32)
 for i,(r,c) in enumerate(pts,1): markers[r,c]=i
 lab=segmentation.watershed(-sm,markers,mask=fg)
 kl=[r.label for r in measure.regionprops(lab) if r.area>=35 and r.bbox[2]<320 and r.bbox[3]<320]
 return len(kl),kl
rows=[]
for _,meta in layout.iterrows():
 im=np.maximum(tifffile.imread(f'/task/input/images/{meta.well}_ch1.tif').astype(float)-D,0)/F
 n,kl=count(im)
 rows.append([meta.well,meta.compound,meta.dose_uM,meta.role,n])
pd.DataFrame(rows,columns=['well','compound','dose','role','n_cells_exact']).to_csv('/task/output/diagnostics/exact_counts.tsv',sep='\t',index=False)
print('controls')
print(pd.DataFrame(rows,columns=['well','compound','dose','role','n']).query("role!='compound'").to_string(index=False))
print('\ncompounds 30->low')
for cp,g in pd.DataFrame(rows,columns=['well','compound','dose','role','n']).query("role=='compound'").groupby('compound'):
 print(f'{cp:10s}',' '.join(f'{r.dose:g}:{r.n}' for _,r in g.sort_values('dose',ascending=False).iterrows()))
