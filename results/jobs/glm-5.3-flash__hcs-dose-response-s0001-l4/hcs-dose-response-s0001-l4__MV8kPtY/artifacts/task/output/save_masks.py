import glob,os,numpy as np,tifffile,multiprocessing as mp
from cellpose import models
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
def work(w):
 x=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32)
 x=np.maximum((x-D)/f,0)
 model=models.Cellpose(gpu=False,model_type='nuclei')
 masks,flows,styles,diam=model.eval(x,diameter=10,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=0)
 np.savez_compressed(f'/task/output/masks_d10_{w}.npz',masks=masks.astype(np.int32),corrected=x,diam=diam)
 return w,int(masks.max()),float(diam)
if __name__=='__main__':
 wells=[os.path.basename(p)[:-4] for p in sorted(glob.glob(base+'/annotations/*.tsv'))]
 with mp.Pool(4) as pool:
  for x in pool.imap_unordered(work,wells): print(x,flush=True)
