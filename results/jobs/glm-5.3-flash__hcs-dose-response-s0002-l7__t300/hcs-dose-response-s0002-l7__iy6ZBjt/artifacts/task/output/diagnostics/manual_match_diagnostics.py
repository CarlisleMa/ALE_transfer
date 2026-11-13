import glob,os,numpy as np,pandas as pd,tifffile
from scipy import ndimage as ndi
from skimage import feature,measure,morphology,segmentation,filters
D=tifffile.imread('/task/input/flatfield/dark.tif').astype(float)
F=tifffile.imread('/task/input/flatfield/flat_ch1.tif').astype(float); F=np.maximum(F-D,0); F/=np.median(F)
rows=[]
for ap in sorted(glob.glob('/task/input/annotations/*.tsv')):
    w=os.path.basename(ap)[:-4]
    im=np.maximum(tifffile.imread(f'/task/input/images/{w}_ch1.tif').astype(float)-D,0)/F
    sm=ndi.gaussian_filter(im,.75)
    pts=feature.peak_local_max(sm,min_distance=3,threshold_abs=900,exclude_border=False)
    pts=pts[~((pts[:,0]>=319)|(pts[:,1]>=319))]
    fg=ndi.binary_fill_holes(morphology.remove_small_objects(sm>filters.threshold_otsu(sm),10))
    markers=np.zeros(im.shape,np.int32)
    for i,(r,c) in enumerate(pts,1): markers[r,c]=i
    lab=segmentation.watershed(-sm,markers,mask=fg)
    regs=[r for r in measure.regionprops(lab) if r.area>=35 and r.bbox[2]<320 and r.bbox[3]<320]
    manual=pd.read_csv(ap,sep='\t'); my=manual.y_px.to_numpy(); mx=manual.x_px.to_numpy()
    assignments=[]
    for r in regs:
        # a manual point matches if inside the region or within 6 px of its centroid
        inside=lab[np.clip(my,0,319).astype(int),np.clip(mx,0,319).astype(int)]==r.label
        near=((my-r.centroid[0])**2+(mx-r.centroid[1])**2)<=36
        nmatch=int(np.sum(inside|near))
        assignments.append((r,nmatch))
    unmatched=[(r,n) for r,n in assignments if n==0]
    matched_multiplicity=pd.Series([n for _,n in assignments if n]).value_counts().sort_index().to_dict()
    print('\n',w,'manual',len(manual),'pipeline',len(regs),'unmatched objects',len(unmatched),'matched object multiplicity',matched_multiplicity)
    print(' unmatched area/minor/major/ecc/bbox:',[(int(r.area),round(float(r.minor_axis_length),1),round(float(r.major_axis_length),1),round(float(r.eccentricity),2),r.bbox) for r,_ in unmatched])
    for r,n in assignments:
        rows.append([w,len(manual),len(regs),r.label,n,int(r.area),float(r.minor_axis_length),float(r.major_axis_length),float(r.eccentricity),r.bbox])
pd.DataFrame(rows,columns=['field','manual','pipeline','label','manual_matches','area','minor','major','ecc','mean_intensity','bbox']).to_csv('/task/output/diagnostics/manual_object_matches.tsv',sep='\t',index=False)
