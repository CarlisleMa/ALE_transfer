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
 inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
 nearest=nuc[inds[0],inds[1]]
 dist=ndi.distance_transform_edt(nuc==0)
 bgs={s:ndi.median_filter(im,size=s,mode='reflect') for s in [31,51,71,91]}
 # ring measurements for each cell
 for i in range(1,int(nuc.max())+1):
  own=nearest==i
  if not own.any(): continue
  cy=ndi.center_of_mass(own)
  yy,xx=np.indices(im.shape)
  rad=np.sqrt((yy-cy[0])**2+(xx-cy[1])**2)
  for r0,r1 in [(1,3),(1,4),(2,4),(2,5),(2,6),(2,8),(4,8),(4,10),(6,10),(6,12),(8,12),(10,15)]:
   sel=own&(rad>=r0)&(rad<r1)&(dist>=1)
   if not sel.any(): continue
   rec=[w,meta.role,i,r0,r1,int(sel.sum())]
   for s,b in bgs.items():
    rec += [float(b[sel].mean()),float(np.clip(im[sel]-b[sel],0,None).mean()),float((im[sel]-b[sel]).mean())]
   rows.append(rec)
cols=['well','role','cell','r0','r1','npx']
for s in [31,51,71,91]: cols += [f'bg{s}',f'clip{s}',f'corr{s}']
d=pd.DataFrame(rows,columns=cols)
d.to_csv('/task/output/diagnostics/local_marker_controls.tsv',sep='\t',index=False)
# aggregate wells
w=d.groupby(['well','role','r0','r1'],as_index=False).mean()
w['row']=w.well.str[0]
print('ROLE MEANS')
print(w.groupby(['r0','r1','role'])[[f'clip{s}' for s in [31,51,71,91]]].mean().to_string(float_format=lambda x:f'{x:.1f}'))
print('\nROW-PAIR REFERENCE/VEHICLE MEDIANS')
for (r0,r1),g in w.groupby(['r0','r1']):
 print(f'\nRING {r0}-{r1}')
 vals=[]
 for _,r in g[g.role=='reference'].iterrows():
  vs=g[(g.role=='vehicle')&(g.row==r.row)]
  z=[r[f'clip{s}']/vs[f'clip{s}'].median() for s in [31,51,71,91]]
  vals.append([r.well,*z])
 print(pd.DataFrame(vals,columns=['well','clip31','clip51','clip71','clip91']).to_string(index=False,float_format=lambda x:f'{x:.2f}'))
