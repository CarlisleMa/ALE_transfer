import numpy as np, pandas as pd, tifffile, os
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
Fl={c:tifffile.imread(f'/task/input/flatfield/flat_ch{c}.tif').astype(float) for c in (1,2)}
for c in Fl: Fl[c]=np.maximum(Fl[c]-D,0)/np.median(Fl[c]-D)
layout=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')
rr,cc=np.indices((320,320))

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
d=pd.read_csv('/task/output/diagnostics/control_radial_raw.tsv',sep='\t')
d=pd.DataFrame(rows,columns=['well','role','compound','dose','cell','r0','r1','n','raw'])
d.to_csv('/task/output/diagnostics/control_radial_raw.tsv',sep='\t',index=False)
# subtract per-cell background based on remote same-territory annuli
for bglo,bghi in [(8,15),(10,20),(12,20),(15,25)]:
 q=d[(d.r0>=bglo)&(d.r1<=bghi)].groupby(['well','cell']).raw.median().rename('bg').reset_index()
 z=d.merge(q,on=['well','cell'],how='left')
 z['corr']=z.raw-z.bg
 print(f'\nBG {bglo}-{bghi}; role means by radial annulus')
 print(z.groupby(['role','r0','r1'])['corr'].mean().unstack([0]).to_string(float_format=lambda x:f'{x:.1f}'))
 print('well control medians for 2-5 and background:')
 s=z[(z.r0==2)&(z.r1==5)].groupby(['well','role']).agg(cy=('corr','mean'),bg=('bg','median')).reset_index()
 print(s.sort_values(['role','well']).to_string(index=False,float_format=lambda x:f'{x:.1f}'))
