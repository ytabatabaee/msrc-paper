"""
This is from here https://github.com/alxsimon/local_pcangsd/blob/main/example.ipynb

"""

import os
import local_pcangsd as lp
import lostruct
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
import seaborn as sns
import dask
import pandas as pd
import numpy as np
import xarray as xr
import sys
from matplotlib.backends.backend_pdf import PdfPages

input = sys.argv[1] 
store = sys.argv[2]

print(input)
print(store)

lp.beagle_to_zarr(input, store, chunksize=100000)

print('finished beagle to zarr!')

ds = lp.load_dataset(store, chunks=100000) # open the Dataset
ds

ds = lp.window(ds, type='position', size=50000, min_variant_number = 500)
# the position argument makes it so that windows are given in their positions
ds

# record the window loci
windf = pd.DataFrame({'window_start' : ds.window_start.to_numpy(),
		      'window_stop' : ds.window_stop.to_numpy(),
		     })

print('starting pca')

#%%time
pca_zarr_store = lp.pca_window(
	ds,
	store=sys.argv[3], # where to store the result
	tmp_folder=sys.argv[4], # need a tmp folder, /tmp/tmp_local_pcangsd is default
	k=10, # number of PCs to retain
)




ds_pca = lp.load_dataset(sys.argv[3])
ds_pca

results = lp.to_lostruct(ds_pca)

print(f"Results on {results.shape[0]} windows")

pc_dists = lostruct.get_pc_dists(results, jax=False)

mds = pcoa(pc_dists)

window_center = lp.get_window_center(ds_pca)

mds.samples['window_center'] = window_center

to_emit = mds.samples[["PC1", "PC2","PC3","PC4","window_center"]]


to_emit.to_csv(sys.argv[8], index=True)

plt.figure()
plt.scatter(x=window_center, y=mds.samples["PC1"])
_ = plt.title("MDS Coordinate 1 (y-axis) compared to Window (x-axis)")
_ = plt.xlabel("Position of window (bp)")
_ = plt.ylabel("MDS 1")

plt.figure()
plt.scatter(x=mds.samples["PC1"], y=mds.samples["PC2"])
_ = plt.xlabel("MDS 1")
_ = plt.ylabel("MDS 2")


propnum=sys.argv[7]
fprop=float(propnum)

mds_12 = mds.samples.loc[:, ["PC1", "PC2"]].copy()
xy = mds_12.to_numpy()
corners = lostruct.corners(xy, prop=fprop) # the proportion here may need to be tuned in some cases. When there are huge inversions, it won't get all of them at 0.05
print(corners[:10])



windf['corner1'] = windf.index.isin(corners[0])
windf['corner2'] = windf.index.isin(corners[1])
windf['corner3'] = windf.index.isin(corners[2])


windf.to_csv(sys.argv[5], index=True)

mds_12['corner'] = 'other'
mds_12['window'] = range(pc_dists.shape[0])
mds_12['window_center'] = window_center
for i in range(3):
	mds_12.iloc[corners[:,i], 2] = f'corner {i+1}'

_ = sns.scatterplot(
	data=mds_12[mds_12.corner=='other'], x="PC1", y="PC2", color='gray',
)
_ = sns.scatterplot(
	data=mds_12[mds_12.corner!='other'], x="PC1", y="PC2", hue="corner",
	hue_order=["corner 1", "corner 2", "corner 3"],
	palette='colorblind',
)



corner1_pca = lp.pcangsd_merged_windows(ds, corners[:,0], k=5)
corner2_pca = lp.pcangsd_merged_windows(ds, corners[:,1], k=5)
corner3_pca = lp.pcangsd_merged_windows(ds, corners[:,2], k=5)

mask1 = ds.windows.isin(corners[:,0])
mask2 = ds.windows.isin(corners[:,1])
mask3 = ds.windows.isin(corners[:,2])


mask = mask1 | mask2 | mask3
mask = ~mask

# use the mask to filter the dataframe
noncorners = ds.windows[mask]


noncorners

baseline_pca = lp.pcangsd_merged_windows(ds, noncorners.values, k=5)




plt.figure()
_ = sns.scatterplot(
	x=corner1_pca[3][0], y=corner1_pca[3][1],)
_ = plt.title("corner 1 pca")
_ = plt.xlabel("PC1")
_ = plt.ylabel("PC2")
_ = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)


plt.figure()
_ = sns.scatterplot(
	x=corner2_pca[3][0], y=corner2_pca[3][1],)
_ = plt.title("corner 2 pca")
_ = plt.xlabel("PC1")
_ = plt.ylabel("PC2")
_ = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

plt.figure()
_ = sns.scatterplot(
	x=corner3_pca[3][0], y=corner3_pca[3][1],)
_ = plt.title("corner 3 pca")
_ = plt.xlabel("PC1")
_ = plt.ylabel("PC2")
_ = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)


plt.figure()
_ = sns.scatterplot(
	x=baseline_pca[3][0], y=baseline_pca[3][1],)
_ = plt.title("baseline pca")
_ = plt.xlabel("PC1")
_ = plt.ylabel("PC2")
_ = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)





pdf = PdfPages(sys.argv[6])
for fig in range(1, plt.gcf().number + 1):
    pdf.savefig( fig, bbox_inches='tight' )
pdf.close()