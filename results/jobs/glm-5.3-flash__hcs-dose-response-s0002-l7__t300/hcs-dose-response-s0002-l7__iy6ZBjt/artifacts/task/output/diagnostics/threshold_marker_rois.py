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
 inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
 nearest=nuc[inds[0],inds[1]]
 rec={'well':w,'role':meta.role,'n_cells':int(nuc.max())}
 for sig in [15,20,30,40]:
  hp=im-ndi.gaussian_filter(im,sig)
  for th in [5,10,15,20,30,40]:
   vals=[]; areas=[]
   for i in range(1,int(nuc.max())+1):
    sel=(nearest==i)&(nuc==0)&(hp>th)
    if sel.any():
     vals.append(float(hp[sel].mean())); areas.append(int(sel.sum()))
    else:
     vals.append(0.0); areas.append(0)
   rec[f'hp{sig}_th{th}']=float(np.mean(vals))
   rec[f'hp{sig}_th{th}_area']=float(np.mean(areas))
 rows.append(rec)
d=pd.DataFrame(rows)
d.to_csv('/task/output/diagnostics/threshold_marker_rois.tsv',sep='\t',index=False)
meas=[c for c in d if c.startswith('hp') and not c.endswith('_area')]
print('ROLE MEDIANS')
print(d.groupby('role')[meas].median().to_string(float_format=lambda x:f'{x:.2f}'))
# spatial normalization row + quadratic column
x=d.well.str[0].map(ord)-65; col=d.well.str[1:].astype(int)-1; d['row']=x;d['col']=col

def design(row,col):
 X=np.zeros((len(row),8))
 for i,r in enumerate(np.asarray(row).astype(int)): X[i,r]=1
 z=(np.asarray(col)-6)/5
 return np.c_[X,z,z*z]
print('\nNORMALIZED REFERENCE RATIOS')
for c in meas:
 v=d[(d.role=='vehicle')&d[c].notna()].reset_index(drop=True)
 p=np.linalg.lstsq(design(v.row,v.col),np.log(np.maximum(v[c].to_numpy(float),1e-6)),rcond=None)[0]
 d[c+'_resp']=d[c]/np.exp(design(d.row,d.col)@p)
 z=d[(d.role=='reference')]
 print(c,' '.join(f'{r.well}:{r[c+"_resp"]:.2f}' for _,r in z.iterrows()))
