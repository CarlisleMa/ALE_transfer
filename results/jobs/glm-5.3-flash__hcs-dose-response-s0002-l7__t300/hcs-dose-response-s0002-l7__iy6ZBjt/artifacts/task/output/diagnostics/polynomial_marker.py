import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
Fl={c:tifffile.imread(f'/task/input/flatfield/flat_ch{c}.tif').astype(float) for c in (1,2)}
for c in Fl: Fl[c]=np.maximum(Fl[c]-D,0)/np.median(Fl[c]-D)
layout=pd.read_csv('/task/input/plate_layout.tsv',sep='\t')
rr,cc=np.indices((320,320))

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

def polyfit_background(im,nuc,degree):
    dist=ndi.distance_transform_edt(nuc==0)
    candidate=(dist>=12)
    # Fit to lower half of candidate pixels, then one robust iteration against positive cell contamination.
    vals=im[candidate]
    cutoff=np.percentile(vals,50)
    use=candidate & (im<=cutoff)
    for _ in range(2):
        z=np.column_stack([ (rr.ravel()/320-.5)**k*(cc.ravel()/320-.5)**j for k in range(degree+1) for j in range(degree+1-k) ])
        coef=np.linalg.lstsq(z[use.ravel()],im.ravel()[use.ravel()],rcond=None)[0]
        bg=(z@coef).reshape(im.shape)
        resid=im-bg
        # keep candidate pixels whose residual is below robust median + 2 MAD
        med=np.median(resid[candidate]); mad=np.median(np.abs(resid[candidate]-med))+1e-6
        use=candidate & (resid<med+2*mad)
    return bg

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
    rec={'well':w,'role':meta.role,'compound':meta.compound,'dose':meta.dose_uM,'n_cells_exact':len(np.unique(nuc))-1}
    for degree in [1,2,3]:
        bg=polyfit_background(im,nuc,degree)
        corrected=im-bg
        for r0,r1 in [(1,4),(2,6),(2,8),(2,10),(3,8),(4,10)]:
            means=[]
            for i in range(1,int(nuc.max())+1):
                sel=(nearest==i)&(nuc==0)&(dist>=r0)&(dist<r1)
                if sel.any(): means.append(float(corrected[sel].mean()))
            rec[f'deg{degree}_{r0}-{r1}']=float(np.mean(means)) if means else np.nan
    rows.append(rec)
d=pd.DataFrame(rows)
d.to_csv('/task/output/diagnostics/polynomial_marker.tsv',sep='\t',index=False)
meas=[c for c in d if c.startswith('deg')]
print('ROLE MEDIANS')
print(d.groupby('role')[meas].median().to_string(float_format=lambda x:f'{x:.2f}'))
# spatial rowcubic vehicle normalization and dose summaries
d['row']=d.well.str[0].map(ord)-65; d['col']=d.well.str[1:].astype(int)-1
qc={'D06','D07','E07','E09','F07'}
def design(row,col):
    r=np.asarray(row)/7; z=(np.asarray(col)-6)/5
    X=np.zeros((len(row),8))
    for i,x in enumerate(np.asarray(row).astype(int)): X[i,x]=1
    return np.c_[X,z,z**2,z**3]
for c in meas:
    v=d[(d.role=='vehicle')&d[c].notna()]
    p=np.linalg.lstsq(design(v.row,v.col),np.log(np.maximum(v[c].to_numpy(float),1e-6)),rcond=None)[0]
    d[c+'_resp']=d[c]/np.exp(design(d.row,d.col)@p)
print('\nREFERENCE NORMALIZED')
for c in meas:
    z=d[(d.role=='reference')&(~d.well.isin(qc))]
    print(c,' '.join(f'{r.well}:{r[c+"_resp"]:.2f}' for _,r in z.iterrows()))
print('\nCOMPOUND NORMALIZED DOSE SERIES')
for c in ['deg2_2-8','deg2_2-10','deg3_2-8','deg3_2-10']:
    print('\n###',c)
    for cp,g in d[(d.role=='compound')&(~d.well.isin(qc))].groupby('compound'):
        h=g.sort_values('dose',ascending=False)
        print(f'{cp:10s}',' '.join(f'{r.dose:g}:{r[c+"_resp"]:.2f}' for _,r in h.iterrows()))
