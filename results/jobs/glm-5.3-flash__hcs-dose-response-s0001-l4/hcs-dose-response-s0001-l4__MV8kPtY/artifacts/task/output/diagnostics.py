import numpy as np, pandas as pd, tifffile, glob, os
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
from skimage.feature import peak_local_max
base='/task/input'
D=tifffile.imread(base+'/flatfield/dark.tif').astype(np.float32)
f1=tifffile.imread(base+'/flatfield/flat_ch1.tif').astype(np.float32)-D; f1/=np.median(f1)
# Full-plate focus metrics on corrected nuclear channel
rows=[]
for p in sorted(glob.glob(base+'/images/*_ch1.tif')):
 w=os.path.basename(p).split('_')[0]
 x=np.maximum((tifffile.imread(p).astype(np.float32)-D)/f1,0)
 # background-normalize using robust lower quantile
 bg=np.quantile(x,.2); y=x/bg
 gy,gx=np.gradient(ndi.gaussian_filter(y,0.5))
 ten=(gx*gx+gy*gy).mean()
 lap=ndi.gaussian_laplace(y,1.0)
 rows.append((w,bg,np.median(x),x.std(),ten,np.abs(lap).mean(),np.quantile(x,.99),np.quantile(x,.999)))
df=pd.DataFrame(rows,columns=['well','bg20','med','std','tenengrad','lapabs','q99','q999']).set_index('well')
df.to_csv('/task/output/focus_metrics.tsv',sep='\t')
for c in df.columns:
 print('\n',c)
 print('min',df[c].idxmin(),df[c].min(),'max',df[c].idxmax(),df[c].max(),'median',df[c].median())
 print('lowest 12:',df[c].nsmallest(12).round(4).to_dict())
 print('highest 5:',df[c].nlargest(5).round(4).to_dict())
