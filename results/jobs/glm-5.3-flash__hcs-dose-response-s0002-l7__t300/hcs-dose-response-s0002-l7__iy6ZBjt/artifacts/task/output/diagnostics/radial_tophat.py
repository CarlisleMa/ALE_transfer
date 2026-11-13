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
 if nuc.max()==0: continue
 # Nearest nucleus territory and distance from all nuclei
 inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
 nearest=nuc[inds[0],inds[1]]
 dist=ndi.distance_transform_edt(nuc==0)
 tophats={f'th{r}':morphology.white_tophat(im,morphology.disk(r)) for r in [4,6,8,10,12,15,20]}
 # Well-level robust background
 bg=np.median(im[dist>14]) if np.any(dist>14) else np.percentile(im,10)
 # Radial bands in territory, raw and top-hat
 for i in range(1,int(nuc.max())+1):
  own=nearest==i
  if not own.any(): continue
  cy=ndi.center_of_mass(own)
  yy,xx=np.indices(im.shape)
  rad=np.sqrt((yy-cy[0])**2+(xx-cy[1])**2)
  for r0,r1 in [(0,2),(2,4),(4,6),(6,8),(8,10),(10,12),(12,15),(15,20),(20,30)]:
   sel=own&(rad>=r0)&(rad<r1)
   if not sel.any(): continue
   rec=[w,meta.role,i,r0,r1,int(sel.sum()),float(im[sel].mean()-bg)]
   rec += [float(v[sel].mean()) for v in tophats.values()]
   rows.append(rec)
cols=['well','role','cell','r0','r1','npx','raw_bgcorr']+[f'th{r}' for r in [4,6,8,10,12,15,20]]
d=pd.DataFrame(rows,columns=cols)
d.to_csv('/task/output/diagnostics/control_radial_tophat.tsv',sep='\t',index=False)
print('saved rows',len(d),'wells',d.well.nunique())
print('\nMEANS BY RADIAL BAND AND ROLE')
print(d.groupby(['r0','r1','role'])[['raw_bgcorr',*[f'th{r}' for r in [4,6,8,10,12,15,20]]]].mean().to_string(float_format=lambda x:f'{x:.1f}'))
print('\nPER-WELL SUMMARY: r=2-6')
s=d[(d.r0==2)&(d.r1==6)].groupby(['well','role'])[['raw_bgcorr',*[f'th{r}' for r in [4,6,8,10,12,15,20]]]].mean().reset_index()
print(s.groupby('role').agg(['mean','median','std']).to_string(float_format=lambda x:f'{x:.1f}'))
print(s.sort_values(['role','well']).to_string(index=False,float_format=lambda x:f'{x:.1f}'))
