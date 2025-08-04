## For GPU-node
salloc --time=02:00:00 -p gpu --gpus-per-node 1 --nodes 1 --ntasks 1 --cpus-per-task 16 --mem 128G

## For CPU-node
salloc -p normal --time=04:00:00 --nodes 1 --ntasks 1 --cpus-per-task 24 --mem 128G

## Setting the directory

cd '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new'
conda activate scrnaseq

## Importing packages
import scanpy as sc
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import numpy as np
import random
import pandas as pd

## Importing cd45 files

cd45=sc.read_h5ad('cd45.h5ad')

## Removing lymphoid cells
cd45=cd45[~cd45.obs['mye_label'].isin(['Lymphoid Lin.'])]

def pseudobulk_sample(adata,sample_key,sample_name,num_replicates1,num_replicates2):
	temp = []
	subset = adata[adata.obs[sample_key].isin([sample_name])]
	if subset.shape[0]>=30:
		subset.X = subset.layers["raw_counts"]
		indices = list(subset.obs_names)
		random.shuffle(indices)
		indices = np.array_split(np.array(indices),num_replicates1)
		for i,pseudo_rep in enumerate(indices):
			rep_adata = sc.AnnData(X=subset[indices[i]].X.sum(axis=0).reshape(1,-1),var=subset[indices[i]].var[[]])
			rep_adata.obs_names = [sample_name + "_" + str(i)]
			temp.append(rep_adata)
		return sc.concat(temp)
	if subset.shape[0]>=20:
		subset.X = subset.layers["raw_counts"]
		indices = list(subset.obs_names)
		random.shuffle(indices)
		indices = np.array_split(np.array(indices),num_replicates2)
		for i,pseudo_rep in enumerate(indices):
			rep_adata = sc.AnnData(X=subset[indices[i]].X.sum(axis=0).reshape(1,-1),var=subset[indices[i]].var[[]])
			rep_adata.obs_names = [sample_name + "_" + str(i)]
			temp.append(rep_adata)
		return sc.concat(temp)
	else:
		subset.X = subset.layers["raw_counts"]
		indices = list(subset.obs_names)
		random.shuffle(indices)
		indices = np.array_split(np.array(indices),1)
		for i,pseudo_rep in enumerate(indices):
			rep_adata = sc.AnnData(X=subset[indices[i]].X.sum(axis=0).reshape(1,-1),var=subset[indices[i]].var[[]])
			rep_adata.obs_names = [sample_name + "_" + str(i)]
			temp.append(rep_adata)
		return sc.concat(temp)

bulk_ls = []

for i in tqdm(cd45.obs["sample"].unique().tolist()):
	bulk_ls.append(pseudobulk_sample(cd45,"sample",i,3,2))

bulk_dat = sc.concat(bulk_ls)
bulk_dat.to_df().to_csv("tissue_pseudobulk.csv")

def mapper(val):
	return val[0:len(val)-2]

bulk_dat.obs['group'] = bulk_dat.obs.index.map(mapper)

metadata = pd.DataFrame(zip(bulk_dat.obs.index,bulk_dat.obs["group"])).rename(columns={0:"sample",1:"group"}).set_index("sample")
metadata.to_csv("tissue_pseudobulk_metadata.csv")

## Setting up deseq2

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

bulk = pd.read_csv("tissue_pseudobulk.csv", index_col=0)
metadata = pd.read_csv("tissue_pseudobulk_metadata.csv",index_col=0)

dds = DeseqDataSet(counts=bulk,metadata=metadata,design_factors="group",refit_cooks=True,n_cpus=6)
sc.pp.filter_genes(dds,min_cells=1)

dds.deseq2()

sc.pp.log1p(dds)
sc.pp.pca(dds,layer='normed_counts')
sc.pl.pca(dds,color=['group'])

## For plotting 3d PCA, let's export the pca coordinates first

pca_df = pd.DataFrame(dds.obsm['X_pca'])
pca_df.index = dds.obs_names
pca_df['group'] = dds.obs['group']
pca_df.to_csv('organ_clustering_PCA_df.csv')

## 3d plot

pca_df = pd.read_csv('organ_clustering_PCA_df.csv',index_col=0)

import distinctipy

colors_ = distinctipy.get_colors(6)
colors_d = dict(zip(pca_df['group'].unique().tolist(),colors_))
pca_df['colors'] = pca_df['group'].map(colors_d)

colors_dd = dict(zip(range(0,18),pca_df['colors']))

## specificying colors now
pca_df=pca_df.reset_index().set_index('group')

for j in tqdm(range(0,360,2)):
	fig = plt.figure(figsize=(5,5))
	ax = plt.axes(projection='3d')
	ax.set_aspect('equal') 
	for i in pca_df.index.unique():
		ax.scatter(pca_df.loc[i]["0"],pca_df.loc[i]["1"],pca_df.loc[i]["2"],c=pca_df.loc[i]['colors'],s=100,label=i)
	xmin = pca_df["0"].min()
	xmax = pca_df["0"].max()
	ymin = pca_df["1"].min()
	ymax = pca_df["1"].max()
	zmin = pca_df["2"].min()
	zmax = pca_df["2"].max()
	ax.plot([xmin,xmax],[ymin,ymin],[zmin,zmin],c='black',lw=1)
	ax.text(xmax,ymin,zmin,'PC1')
	ax.plot([xmin,xmin],[ymin,ymax],[zmin,zmin],c='black',lw=1)
	ax.text(xmin,ymax,zmin,'PC2')
	ax.plot([xmin,xmin],[ymin,ymin],[zmin,zmax],c='black',lw=1)
	ax.text(xmin,ymin,zmax,'PC3')
	ax.view_init(20,j)
	plt.legend()
	ax.axis('off')
	plt.tight_layout()
	plt.savefig('organ_3d_views/pngs/3d_pca_'+str(j)+'.png')
	plt.savefig('organ_3d_views/svgs/3d_pca_'+str(j)+'.svg')
	plt.clf()
	plt.close('all')


