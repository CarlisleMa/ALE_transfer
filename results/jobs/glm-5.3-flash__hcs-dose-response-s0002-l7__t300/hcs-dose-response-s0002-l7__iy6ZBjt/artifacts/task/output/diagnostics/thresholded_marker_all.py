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

def background(im,nuc):
    dist=ndi.distance_transform_edt(nuc==0)
    candidate=dist>=12
    use=candidate&(im<=np.percentile(im[candidate],50))
    z=np.column_stack([(rr.ravel()/320-.5)**k*(cc.ravel()/320-.5)**j for k in range(3) for j in range(3-k)])
    for _ in range(2):
        coef=np.linalg.lstsq(z[use.ravel()],im.ravel()[use.ravel()],rcond=None)[0]
        bg=(z@coef).reshape(im.shape)
        residual=im-bg
        med=np.median(residual[candidate]); mad=np.median(np.abs(residual[candidate]-med))+1e-6
        use=candidate&(residual<med+2*mad)
    return bg,med,mad

rows=[]
for _,meta in layout.iterrows():
    w=meta.well
    nuc=nuclear_segment(np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/Fl[1])
    im=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch2.tif').astype(float)-D,0)/Fl[2]
    if nuc.max()==0: continue
    bg,med,mad=background(im,nuc)
    corrected=im-bg
    inds=ndi.distance_transform_edt(nuc==0,return_indices=True)[1]
    nearest=nuc[inds[0],inds[1]]
    dist=ndi.distance_transform_edt(nuc==0)
    rec={'well':w,'role':meta.role,'compound':meta.compound,'dose':meta.dose_uM,'n_cells_exact':len(np.unique(nuc))-1,'bg_med':med,'bg_mad':mad}
    thresholds={'mad3':max(10,med+3*1.4826*mad),'fixed15':15,'fixed30':30,'fixed60':60}
    for name,th in thresholds.items():
        means={}
        for maxdist in [12,15,20,30]:
            vals=[]; areas=[]
            for i in range(1,int(nuc.max())+1):
                sel=(nearest==i)&(nuc==0)&(dist<=maxdist)&(corrected>th)
                if sel.any():
                    vals.append(float(corrected[sel].mean())); areas.append(int(sel.sum()))
                else:
                    vals.append(0.0); areas.append(0)
            means[f'th_{name}_d{maxdist}']=float(np.mean(vals))
            means[f'th_{name}_d{maxdist}_area']=float(np.mean(areas))
        rec.update(means)
    rows.append(rec)
d=pd.DataFrame(rows)
d.to_csv('/task/output/diagnostics/thresholded_marker_all.tsv',sep='\t',index=False)
cols=[c for c in d if c.startswith('th_') and not c.endswith('_area')]
print('ROLE MEDIANS')
print(d.groupby('role')[cols].median().to_string(float_format=lambda x:f'{x:.2f}'))
# row-specific linear/quadratic and global row+quadratic normalizations
d['row']=d.well.str[0].map(ord)-65; d['col']=d.well.str[1:].astype(int)-1
qc={'D06','D07','E07','E09','F07'}
for c in cols:
    # global row fixed effects + quadratic column on log values
    v=d[(d.role=='vehicle')&d[c].notna()]
    X=np.zeros((len(v),8))
    for i,r in enumerate(v.row.astype(int)): X[i,r]=1
    X=np.c_[X,(v.col-6)/5,((v.col-6)/5)**2]
    p=np.linalg.lstsq(X,np.log(np.maximum(v[c].to_numpy(float),1e-6)),rcond=None)[0]
    Xa=np.zeros((len(d),8))
    for i,r in enumerate(d.row.astype(int)): Xa[i,r]=1
    Xa=np.c_[Xa,(d.col-6)/5,((d.col-6)/5)**2]
    d[c+'_resp']=d[c]/np.exp(Xa@p)
print('\nREFERENCE NORMALIZED')
for c in cols:
    z=d[(d.role=='reference')&(~d.well.isin(qc))]
    print(c,' '.join(f'{r.well}:{r[c+"_resp"]:.2f}' for _,r in z.iterrows()))
print('\nSELECTED COMPOUND SERIES')
for c in ['th_mad3_d15','th_fixed15_d15','th_fixed30_d15','th_fixed60_d15','th_fixed30_d20']:
    print('\n###',c)
    for cp,g in d[(d.role=='compound')&(~d.well.isin(qc))].groupby('compound'):
        h=g.sort_values('dose',ascending=False)
        print(f'{cp:10s}',' '.join(f'{r.dose:g}:{r[c+"_resp"]:.2f}' for _,r in h.iterrows()))
