import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
Fl={c:tifffile.imread(f'/task/input/flatfield/flat_ch{c}.tif').astype(float) for c in (1,2)}
for c in Fl: Fl[c]=np.maximum(Fl[c]-D,0)/np.median(Fl[c]-D)
layout=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')

def nuclear_segment(im):
 sm=ndi.gaussian_filter(im,0.75)
 pts=feature.peak_local_max(sm,min_distance=3,threshold_abs=900,exclude_border=False)
 pts=pts[~((pts[:,0]>=319)|(pts[:,1]>=319))]
 fg=ndi.binary_fill_holes(morphology.remove_small_objects(sm>filters.threshold_otsu(sm),10))
 markers=np.zeros(im.shape,np.int32)
 for i,(r,c) in enumerate(pts,1): markers[r,c]=i
 lab=segmentation.watershed(-sm,markers,mask=fg)
 kl=[r.label for r in measure.regionprops(lab) if r.area>=35 and r.bbox[2]<320 and r.bbox[3]<320]
 return np.where(np.isin(lab,kl),lab,0)

rows=[]
for _,meta in layout[layout.role.isin(['vehicle','reference'])].iterrows():
 w=meta.well
 nuc=nuclear_segment(np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/Fl[1])
 im=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch2.tif').astype(float)-D,0)/Fl[2]
 rec={'well':w,'role':meta.role,'n_cells':int(nuc.max())}
 for sig in [10,15,20,30]:
  bg=ndi.gaussian_filter(im,sig)
  hp=im-bg
  ot=filters.threshold_otsu(hp)
  fg=ndi.binary_fill_holes(morphology.remove_small_objects(hp>ot,20))
  # Include nuclei in mask so every counted cell has an ROI; watershed nuclear seeds.
  mask=fg | (nuc>0)
  lab=segmentation.watershed(-hp,nuc,mask=mask)
  means=[]; areas=[]
  for i in range(1,int(nuc.max())+1):
   sel=lab==i
   if sel.any():
    means.append(float(np.clip(hp[sel],0,None).mean()))
    areas.append(int(sel.sum()))
  rec[f'obj_sig{sig}']=np.mean(means) if means else np.nan
  rec[f'obj_sig{sig}_area']=np.mean(areas) if areas else np.nan
  rec[f'ot_sig{sig}']=ot
  rec[f'fg_sig{sig}']=fg.mean()
 rows.append(rec)
d=pd.DataFrame(rows)
d.to_csv('/task/output/diagnostics/marker_object_rois.tsv',sep='\t',index=False)
print(d[['well','role','n_cells','fg_sig10','fg_sig15','fg_sig20','fg_sig30','obj_sig10','obj_sig15','obj_sig20','obj_sig30']].to_string(index=False,float_format=lambda x:f'{x:.2f}'))
print('\nROLE SUMMARY')
print(d.groupby('role')[['obj_sig10','obj_sig15','obj_sig20','obj_sig30']].agg(['mean','median','std']).to_string(float_format=lambda x:f'{x:.2f}'))
