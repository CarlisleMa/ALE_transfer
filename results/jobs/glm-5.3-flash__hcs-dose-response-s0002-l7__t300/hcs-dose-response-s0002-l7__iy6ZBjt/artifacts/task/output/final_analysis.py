import os, glob, math
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from skimage import feature, measure, morphology, segmentation, filters

IN='/task/input'
OUT='/task/output'
os.makedirs(OUT,exist_ok=True)

# Dark/flat-field correction
dark=tifffile.imread(f'{IN}/flatfield/dark.tif').astype(np.float32)
flat={}
for ch in (1,2):
    f=tifffile.imread(f'{IN}/flatfield/flat_ch{ch}.tif').astype(np.float32)
    f=np.maximum(f-dark,0)
    flat[ch]=f/np.median(f)

layout=pd.read_csv(f'{IN}/plate_layout.tsv',sep='\t')
rr,cc=np.indices((320,320))

# Validated nuclear segmentation/counting procedure.
# Intact nuclei are enforced by a 35-pixel minimum area; bottom/right-cut objects are excluded.
def nuclear_segment(im):
    sm=ndi.gaussian_filter(im,0.75)
    pts=feature.peak_local_max(sm,min_distance=3,threshold_abs=900,exclude_border=False)
    pts=pts[~((pts[:,0]>=319)|(pts[:,1]>=319))]
    fg=ndi.binary_fill_holes(morphology.remove_small_objects(sm>filters.threshold_otsu(sm),10))
    markers=np.zeros(im.shape,dtype=np.int32)
    for i,(r,c) in enumerate(pts,1):
        markers[r,c]=i
    lab=segmentation.watershed(-sm,markers,mask=fg)
    kept=[]
    for r in measure.regionprops(lab):
        minr,minc,maxr,maxc=r.bbox
        if r.area>=35 and maxr<320 and maxc<320:
            kept.append(r.label)
    nuc=np.where(np.isin(lab,kept),lab,0).astype(np.int32)
    return nuc,len(kept)

# Focus score is measured on nuclear edges, so it is not merely a cell-density metric.
def focus_score(im):
    sm=ndi.gaussian_filter(im,0.75)
    fg=sm>filters.threshold_otsu(sm)
    zone=ndi.binary_dilation(fg) & ~ndi.binary_erosion(fg)
    gy,gx=np.gradient(ndi.gaussian_filter(im,0.75))
    edge=float(np.mean(np.hypot(gx,gy)[zone]))
    return edge/np.sqrt(max(float(fg.mean()),1e-12))

# Robust quadratic background fitted to pixels distant from nuclei.
def marker_background(im,nuc):
    dist=ndi.distance_transform_edt(nuc==0)
    candidate=dist>=12
    if not candidate.any():
        return np.full_like(im,np.percentile(im,10))
    use=candidate & (im<=np.percentile(im[candidate],50))
    z=np.column_stack([
        (rr.ravel()/320.0-0.5)**i*(cc.ravel()/320.0-0.5)**j
        for i in range(3) for j in range(3-i)
    ])
    y=im.ravel()
    med=float(np.median(y[candidate.ravel()]))
    mad=1.0
    for _ in range(2):
        coef=np.linalg.lstsq(z[use.ravel()],y[use.ravel()],rcond=None)[0]
        bg=(z@coef).reshape(im.shape)
        resid=y-bg.ravel()
        med=float(np.median(resid[candidate.ravel()]))
        mad=float(np.median(np.abs(resid[candidate.ravel()]-med)))+1e-6
        use=candidate & (resid.reshape(im.shape)<med+2.0*mad)
    return bg

rows=[]
for _,meta in layout.iterrows():
    w=meta['well']
    im1=np.maximum(tifffile.imread(f'{IN}/images/{w}_ch1.tif').astype(np.float32)-dark,0)/flat[1]
    im2=np.maximum(tifffile.imread(f'{IN}/images/{w}_ch2.tif').astype(np.float32)-dark,0)/flat[2]
    nuc,n_cells=nuclear_segment(im1)
    score=focus_score(im1)
    # Sharp focused fields are far above 300; the five visibly blurred fields are far below.
    qc_pass=int(score>=300.0)

    marker=np.nan
    marker_cells=[]
    if n_cells>0:
        bg=marker_background(im2,nuc)
        corrected=im2-bg
        _,indices=ndi.distance_transform_edt(nuc==0,return_indices=True)
        nearest=nuc[indices[0],indices[1]]
        dist=ndi.distance_transform_edt(nuc==0)
        for label in range(1,int(nuc.max())+1):
            if not np.any(nuc==label):
                continue
            sel=(nearest==label)&(nuc==0)&(dist>=2)&(dist<8)
            if sel.any():
                marker_cells.append(float(corrected[sel].mean()))
        if marker_cells:
            marker=float(np.mean(marker_cells))

    rows.append(dict(
        well=w, compound=meta['compound'], dose_uM=meta['dose_uM'], role=meta['role'],
        n_cells=int(n_cells), marker_per_cell=marker, qc_pass=qc_pass,
        focus_score=score, n_marker_cells=len(marker_cells)
    ))

df=pd.DataFrame(rows)
df['row']=df.well.str[0].map(ord)-65
df['col']=df.well.str[1:].astype(int)-1

# Same-position vehicle-only counterfactuals from the distributed vehicle controls.
def design(row,col):
    row=np.asarray(row,dtype=int); col=np.asarray(col,dtype=float)
    X=np.zeros((len(row),8),dtype=float)
    for i,r in enumerate(row):
        X[i,r]=1.0
    z=(col-6.0)/5.0
    return np.c_[X,z,z*z]

for target,log_scale in [('n_cells',False),('marker_per_cell',True)]:
    veh=df[(df.role=='vehicle') & df[target].notna()].copy()
    y=veh[target].to_numpy(float)
    fit_y=np.log(np.maximum(y,1e-12)) if log_scale else y
    coef=np.linalg.lstsq(design(veh.row,veh.col),fit_y,rcond=None)[0]
    pred=design(df.row,df.col)@coef
    if log_scale:
        pred=np.exp(pred)
    df[f'vehicle_pred_{target}']=pred
    df[f'response_{target}']=df[target]/pred

# SOP four-parameter logistic: fixed upper asymptote at one, 0<=bottom<1, h>0.
# Global optimization avoids the narrow-slope local minima seen with local starting points.
from scipy.optimize import differential_evolution

def logistic(c,bottom,ec50,h):
    c=np.asarray(c,float)
    z=np.clip(h*np.log(np.maximum(c,1e-300)/max(ec50,1e-300)),-50,50)
    return bottom+(1.0-bottom)/(1.0+np.exp(z))

def fit_response(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    def objective(p):
        bottom,log_ec50,h=p
        return float(np.sum((logistic(x,bottom,np.exp(log_ec50),h)-y)**2))
    res=differential_evolution(
        objective,
        bounds=[(0.0,0.999999),(np.log(1e-12),np.log(1e6)),(0.05,100.0)],
        seed=12345,tol=1e-9,popsize=20,maxiter=500,polish=True
    )
    bottom=float(res.x[0]); ec50=float(np.exp(res.x[1])); h=float(res.x[2])
    if bottom<0.5:
        potency=float(ec50*(1.0/(1.0-2.0*bottom))**(1.0/h))
    else:
        potency=np.nan
    return dict(bottom=bottom,ec50=ec50,h=h,potency=potency,sse=float(res.fun))

usable=df[df.qc_pass==1]
potencies=[]
for compound,g in usable[usable.role=='compound'].groupby('compound',sort=False):
    g=g.sort_values('dose_uM')
    marker_fit=fit_response(g.dose_uM,g.response_marker_per_cell)
    count_fit=fit_response(g.dose_uM,g.response_n_cells)
    def result(fit):
        if fit is None:
            return 'inactive',fit
        pot=fit['potency']
        # SOP: a fitted 0.5 crossing outside the tested range is inactive.
        if not np.isfinite(pot) or pot<g.dose_uM.min() or pot>g.dose_uM.max():
            return 'inactive',fit
        return pot,fit
    ic50,mfit=result(marker_fit)
    gi50,cfit=result(count_fit)
    potencies.append(dict(compound=compound,ic50_marker_uM=ic50,gi50_count_uM=gi50,
                          marker_fit=marker_fit,count_fit=count_fit))
pot=pd.DataFrame(potencies)

# Required deliverables.
well_out=df[['well','n_cells','marker_per_cell','qc_pass','focus_score',
             'vehicle_pred_n_cells','response_n_cells',
             'vehicle_pred_marker_per_cell','response_marker_per_cell']].copy()
well_out['n_cells']=well_out.n_cells.astype(int)
well_out['qc_pass']=well_out.qc_pass.astype(int)
well_out.to_csv(f'{OUT}/well_measurements.tsv',sep='\t',index=False,na_rep='NA',float_format='%.6f')

validation=[]
for _,meta in layout[layout.role.isin(['compound','vehicle'])].iterrows():
    ap=f'{IN}/annotations/{meta.well}.tsv'
    if os.path.exists(ap):
        manual=int(len(pd.read_csv(ap,sep='\t')))
        pipeline=int(df.loc[df.well==meta.well,'n_cells'].iloc[0])
        validation.append(dict(field=meta.well,manual_count=manual,pipeline_count=pipeline))
pd.DataFrame(validation,columns=['field','manual_count','pipeline_count']).to_csv(f'{OUT}/validation.tsv',sep='\t',index=False)

# Main required columns, with fit diagnostics in extra columns.
dose_out=pot[['compound','ic50_marker_uM','gi50_count_uM']].copy()
for rd in ('marker','count'):
    for key in ('bottom','ec50','h','potency','sse'):
        vals=[]
        for _,r in pot.iterrows():
            fit=r[f'{rd}_fit']
            vals.append(np.nan if fit is None else fit.get(key,np.nan))
        dose_out[f'{rd}_{key}']=vals
dose_out.to_csv(f'{OUT}/dose_response.tsv',sep='\t',index=False,na_rep='NA',float_format='%.8g')

# Human-readable processing summary (not a required deliverable).
print('Wells:',len(df),'QC failures:',df.loc[df.qc_pass==0,'well'].tolist())
print('\nValidation:')
print(pd.DataFrame(validation).to_string(index=False))
print('\nPotencies:')
print(dose_out[['compound','ic50_marker_uM','gi50_count_uM']].to_string(index=False))
