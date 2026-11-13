import numpy as np,pandas as pd,tifffile,glob,os
from scipy import ndimage as ndi
from skimage.measure import regionprops
d=pd.read_csv('/task/output/final_object_diagnostics.tsv',sep='\t')
# Final counting candidates will be applied later; prototype using masks and final nuclei
for p in sorted(glob.glob('/task/input/annotations/*.tsv')):
 w=os.path.basename(p)[:-4]; z=np.load(f'/task/output/final_masks_{w}.npz'); masks=z['masks']; ch2=z['ch2']; ch1=z['ch1']
 # SOP counting frame and validated fragment rule (area>=50)
 nuc=np.zeros_like(masks,dtype=bool); keep=[]
 for r in regionprops(masks):
  if r.bbox[2]<320 and r.bbox[3]<320 and r.area>=50:
   nuc|=masks==r.label; keep.append(r.label)
 # distance to any counted nucleus and nearest-label expansion
 dist,inds=ndi.distance_transform_edt(~nuc,return_indices=True)
 nearest=np.zeros_like(masks)
 nearest[nuc]=masks[nuc]
 bgmask=~nuc
 nearest[bgmask]=masks[inds[0][bgmask],inds[1][bgmask]]
 # global robust background from low marker pixels
 hist,edges=np.histogram(ch2,bins=256,range=(0,1024)); bgbin=np.argmax(hist); bgmode=(edges[bgbin]+edges[bgbin+1])/2
 profiles=[]
 for lab in keep:
  r=next(r for r in regionprops(masks,intensity_image=ch1) if r.label==lab)
  cy,cx=r.centroid
  assign=(nearest==lab)&(~nuc)
  for rad in (8,10,12,14,16):
   territory=assign&(dist<=rad)
   profiles.append((w,lab,rad,territory.sum(),ch2[territory].mean(),ch2[territory].mean()-bgmode,np.quantile(ch2[territory],.1),np.quantile(ch2[territory],.5)))
 prof=pd.DataFrame(profiles,columns=['well','label','radius','npx','mean','bgsub','q10','med'])
 print('\n',w,'global bg mode',round(bgmode,1),'cells',len(keep))
 print(prof.groupby('radius')[['npx','mean','bgsub','q10','med']].median().round(1).to_string())
