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
for _,meta in layout.iterrows():
 w=meta.well
 nuc=nuclear_segment(np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/Fl[1])
 im=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch2.tif').astype(float)-D,0)/Fl[2]
 if nuc.max()==0:
  continue
 inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
 nearest=nuc[inds[0],inds[1]]
 dist=ndi.distance_transform_edt(nuc==0)
 # fast local background surfaces
 bgs={f'g{s}':ndi.gaussian_filter(im,s) for s in [2,3,5,8,12,20]}
 # per-cell ring means for raw-minus-global and highpass images
 rec={'well':w,'role':meta.role,'compound':meta.compound,'dose':meta.dose_uM,'n_cells':int(nuc.max())}
 for r0,r1 in [(1,3),(1,4),(2,4),(2,5),(2,6),(2,8),(3,6),(4,8)]:
  means={k:[] for k in ['raw'] + list(bgs)}
  for i in range(1,int(nuc.max())+1):
   sel=(nearest==i)&(dist>=r0)&(dist<r1)
   if not sel.any(): continue
   means['raw'].append(float(im[sel].mean()))
   for k,b in bgs.items(): means[k].append(float(np.clip(im[sel]-b[sel],0,None).mean()))
  for k,v in means.items():
   rec[f'{k}_{r0}-{r1}']=float(np.mean(v)) if v else np.nan
 # remote extracellular background per cell: 7-12 px from nucleus, own Voronoi territory
 remote=[]; remote_med=[]
 for i in range(1,int(nuc.max())+1):
  sel=(nearest==i)&(dist>=7)&(dist<=12)
  if sel.any(): remote.append(float(im[sel].mean())); remote_med.append(float(np.median(im[sel])))
 rec['bg_remote_mean']=np.mean(remote) if remote else np.nan
 rec['bg_remote_median']=np.median(remote_med) if remote_med else np.nan
 for r0,r1 in [(1,3),(1,4),(2,4),(2,5),(2,6),(2,8)]:
  means=[]
  for i in range(1,int(nuc.max())+1):
   cy=(nearest==i)&(dist>=r0)&(dist<r1)
   bg=(nearest==i)&(dist>=7)&(dist<=12)
   if cy.any() and bg.any(): means.append(float(im[cy].mean()-np.median(im[bg])))
  rec[f'remote_{r0}-{r1}']=np.mean(means) if means else np.nan
 rows.append(rec)
pd.DataFrame(rows).to_csv('/task/output/diagnostics/fast_marker_candidates.tsv',sep='\t',index=False)
print('saved',len(rows))
