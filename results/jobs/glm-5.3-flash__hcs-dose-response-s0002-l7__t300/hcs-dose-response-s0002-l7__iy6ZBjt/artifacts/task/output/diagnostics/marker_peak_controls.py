import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
Fl={c:tifffile.imread(f'/task/input/flatfield/flat_ch{c}.tif').astype(float) for c in (1,2)}
for c in Fl: Fl[c]=np.maximum(Fl[c]-D,0)/np.median(Fl[c]-D)
layout=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')

def nuclear_segment(im):
    sm=ndi.gaussian_filter(im,.75)
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
    dist=ndi.distance_transform_edt(nuc==0)
    # local smooth background and small-structure top-hats
    bg=ndi.gaussian_filter(im,20)
    hp=im-bg
    th={r:morphology.white_tophat(im,morphology.disk(r)) for r in [2,3,4,6,8]}
    rec={'well':w,'role':meta.role,'n_cells':int(nuc.max())}
    for radius in [4,6,8,10,12]:
        vals=[]
        for i in range(1,int(nuc.max())+1):
            sel=(nearest==i)&(nuc==0)&(dist<=radius)
            if sel.any():
                a=hp[sel]
                vals.append([np.mean(a),np.percentile(a,75),np.percentile(a,90),np.max(a)]+[float(v[sel].mean()) for v in th.values()])
        if vals:
            arr=np.asarray(vals)
            for j,nm in enumerate(['mean','p75','p90','max']): rec[f'hp_{nm}_r{radius}']=float(np.mean(arr[:,j]))
            for j,r in enumerate([2,3,4,6,8]): rec[f'th{r}_r{radius}']=float(np.mean(arr[:,4+j]))
    rows.append(rec)
d=pd.DataFrame(rows)
d.to_csv('/task/output/diagnostics/marker_peak_controls.tsv',sep='\t',index=False)
print(d.groupby('role')[[c for c in d if c.startswith(('hp_','th'))]].median().to_string(float_format=lambda x:f'{x:.2f}'))
# Normalize within row to vehicle median and print reference ratios
v=d[d.role=='vehicle']; print('\nREFERENCE / ROW VEHICLE MEDIAN')
for _,r in d[d.role=='reference'].iterrows():
    q=v[v.well.str[0]==r.well[0]]
    print(r.well,' '.join(f'{c}:{r[c]/q[c].median():.2f}' for c in [x for x in d if x.startswith('hp_mean')]+[x for x in d if x.startswith('th') and x.endswith('r8')]))
