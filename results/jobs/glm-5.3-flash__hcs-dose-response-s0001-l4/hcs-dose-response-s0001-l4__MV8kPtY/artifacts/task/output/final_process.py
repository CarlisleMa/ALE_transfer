import os, numpy as np, pandas as pd, tifffile, multiprocessing as mp
from scipy import ndimage as ndi
from skimage.measure import regionprops
from cellpose import models
import torch

BASE='/task/input'
DARK=tifffile.imread(BASE+'/flatfield/dark.tif').astype(np.float32)
FL=[]
for ch in (1,2):
    f=tifffile.imread(f'{BASE}/flatfield/flat_ch{ch}.tif').astype(np.float32)-DARK
    f/=np.median(f)
    FL.append(f)
NC=CC=None

def init_models():
    global NC,CC
    torch.set_num_threads(1)
    NC=models.CellposeModel(gpu=False,model_type='nuclei')
    CC=models.CellposeModel(gpu=False,model_type='cyto3')

def correct(w):
    out=[]
    for ch in (1,2):
        raw=tifffile.imread(f'{BASE}/images/{w}_ch{ch}.tif').astype(np.float32)
        out.append(np.maximum((raw-DARK)/FL[ch-1],0))
    return out

def process(well):
    ch1,ch2=correct(well)
    # Final validated nuclear counting settings.
    nm,_,_=NC.eval(ch1,diameter=8.5,channels=[0,0],normalize=True,
                   cellprob_threshold=-1,flow_threshold=2,min_size=15)
    counted=[]
    for r in regionprops(nm):
        # SOP counting frame: top/left cuts are kept; bottom/right cuts are dropped.
        if r.bbox[2] < 320 and r.bbox[3] < 320 and r.area >= 50:
            counted.append(int(r.label))
    n_cells=len(counted)
    counted_mask=np.isin(nm,counted)
    nuclear_any=nm>0

    # Final cytoplasmic territories from cyto3, with the marker as cytoplasm and ch1 as nucleus.
    cm,_,_=CC.eval(np.stack([ch1,ch2]),channels=[2,1],diameter=15,normalize=True,
                   flow_threshold=0.6,cellprob_threshold=0,min_size=30)
    cell_foreground=cm>0

    # Assign every extracellular/cytoplasmic pixel to its nearest counted nucleus.
    if n_cells:
        dist_to_counted,nearest_indices=ndi.distance_transform_edt(~counted_mask,return_indices=True)
        nearest_counted=np.where(counted_mask,nm,nm[nearest_indices[0],nearest_indices[1]])
    else:
        dist_to_counted=None
        nearest_counted=np.zeros_like(nm)

    # Extracellular background is estimated outside both cyto3 cells and nuclear objects.
    background_mask=(cm==0)&(~nuclear_any)
    if background_mask.any():
        background=float(np.median(ch2[background_mask]))
    else:
        background=float(np.percentile(ch2,10))

    per_cell=[]
    fallback=0
    if n_cells:
        # Each counted nucleus takes all cyto3 cells it overlaps, excluding nuclear pixels.
        # Nearest-counted-nucleus partitioning handles cyto3 masks containing multiple nuclei.
        for label in counted:
            nucleus=nm==label
            overlap_labels=np.unique(cm[nucleus])
            overlap_labels=overlap_labels[overlap_labels>0]
            cyto=(nearest_counted==label)&(~nuclear_any)
            if len(overlap_labels):
                cyto &= np.isin(cm,overlap_labels)
            if not cyto.any():
                # Local fallback for a rare nucleus with no cyto3-derived cytoplasmic pixels.
                cyto=(nearest_counted==label)&(~nuclear_any)&(dist_to_counted<=10)
                fallback += 1
            if cyto.any():
                per_cell.append(float(ch2[cyto].mean()-background))
            else:
                per_cell.append(0.0)

    # Focus metric on the corrected nuclear channel. The five low-focus fields are
    # separated from all other fields by a large gap in this metric.
    gy,gx=np.gradient(ndi.gaussian_filter(ch1,0.5))
    focus=float(np.mean(gx*gx+gy*gy))
    qc=int(focus >= 2000.0)
    marker=float(np.mean(per_cell)) if per_cell else np.nan
    return dict(well=well,n_cells=n_cells,marker_per_cell_raw=marker,
                marker_per_cell_bgsub=marker,
                background=background,n_cyto_fallback=fallback,
                focus_ten=focus,qc_pass=qc)

def main():
    layout=pd.read_csv(BASE+'/plate_layout.tsv',sep='\t')
    rows=[]
    with mp.Pool(4,initializer=init_models) as pool:
        for i,row in enumerate(pool.imap_unordered(process,layout.well.tolist()),1):
            rows.append(row)
            print('processed',i,row['well'],flush=True)
    out=pd.DataFrame(rows).merge(layout,on='well',how='left').sort_values('well')
    out.to_csv('/task/output/final_measurements_raw.tsv',sep='\t',index=False)
    print('\nQC failures:',out.loc[out.qc_pass==0,'well'].tolist())
    print('count range:',out.n_cells.min(),out.n_cells.max())
    print('marker NA wells:',out.loc[out.marker_per_cell_bgsub.isna(),'well'].tolist())

if __name__=='__main__':
    main()
