import numpy as np,pandas as pd,tifffile,time,multiprocessing as mp
from cellpose import models
from skimage.measure import regionprops
base='/task/input'; D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f/=np.median(f)
wells=['A06','A11','A12','C12','D10','E05','G07','H12']
def work(a):
 diam,cp,w=a
 x=tifffile.imread(f'{base}/images/{w}_ch1.tif').astype(np.float32); x=np.maximum((x-D)/f,0)
 model=models.Cellpose(gpu=False,model_type='nuclei')
 masks,_,_,_=model.eval(x,diameter=diam,channels=[0,0],normalize=True,flow_threshold=.6,cellprob_threshold=cp)
 regs=regionprops(masks)
 # edge frame rule: exclude max row/col at 320
 count=sum(r.bbox[2]<320 and r.bbox[3]<320 for r in regs)
 areas=np.array([r.area for r in regs if r.bbox[2]<320 and r.bbox[3]<320])
 manual=len(pd.read_csv(f'{base}/annotations/{w}.tsv',sep='\t'))
 return diam,cp,w,count,manual,int((areas>=30).sum()),int((areas>=40).sum()),np.percentile(areas,[0,5,10,25,50,90,100]).round(1).tolist()
if __name__=='__main__':
 jobs=[]
 for diam in (8,8.5,9,10,11):
  for cp in (-1,-0.5,0,0.5):
   for w in wells: jobs.append((diam,cp,w))
 with mp.Pool(4) as pool:
  for x in pool.imap_unordered(work,jobs): print(x,flush=True)
