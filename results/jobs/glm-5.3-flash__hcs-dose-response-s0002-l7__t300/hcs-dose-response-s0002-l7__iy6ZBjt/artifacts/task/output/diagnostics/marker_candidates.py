import os,numpy as np,pandas as pd,tifffile
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
for _,meta in layout.iterrows():
 w=meta.well
 nuc=nuclear_segment(np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/Fl[1])
 im2=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch2.tif').astype(float)-D,0)/Fl[2]
 if nuc.max()==0:
  rows.append([w,meta.compound,meta.dose_uM,meta.role,0,*[np.nan]*12]); continue
 # nearest nucleus territory
 inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
 nearest=np.zeros_like(nuc); nearest[nuc>0]=nuc[nuc>0]; nearest=nearest[inds[0],inds[1]]
 nucdist=ndi.distance_transform_edt(nuc==0)
 # Background estimates: robust low values and smoothed morphological opening
 bgs={
  'p10':np.percentile(im2,10),
  'p20':np.percentile(im2,20),
  'open':np.median(ndi.gaussian_filter(ndi.grey_opening(im2,size=(31,31)),3)),
  'far':np.median(im2[nucdist>12]) if np.any(nucdist>12) else np.percentile(im2,10),
 }
 # Marker foreground: threshold background + margin, moderate smoothing
 vals={}
 for margin in [15,25,40]:
  th=bgs['open']+margin
  fg2=ndi.binary_fill_holes(ndi.gaussian_filter(im2,1)>th)
  assigned=(fg2)&(nucdist>=1)&(nucdist<=18)&(nearest>0)
  # cap remote pixels by Voronoi assignment but exclude ambiguous pixels farther than 18px
  means=[]
  for i in range(1,int(nuc.max())+1):
   sel=assigned&(nearest==i)
   if sel.any(): means.append(np.mean(im2[sel])-bgs['open'])
  vals[f'mk_open{margin}']=np.mean(means) if means else np.nan
  vals[f'npx_open{margin}']=len(means)
 # ring candidates
 for outer in [3,4,5,6,8]:
  ring=(nucdist>=2)&(nucdist<=outer)&(nearest>0)
  means=[]
  for i in range(1,int(nuc.max())+1):
   sel=ring&(nearest==i)
   if sel.any(): means.append(np.mean(im2[sel])-bgs['open'])
  vals[f'ring{outer}']=np.mean(means) if means else np.nan
 rows.append([w,meta.compound,meta.dose_uM,meta.role,len(np.unique(nuc))-1,bgs['p10'],bgs['p20'],bgs['open'],bgs['far'],*list(vals.values())])
cols=['well','compound','dose','role','n_cells','bg_p10','bg_p20','bg_open','bg_far']+[x for x in vals]
pd.DataFrame(rows,columns=cols).to_csv('/task/output/diagnostics/marker_candidates.tsv',sep='\t',index=False)
print('saved',len(rows))
