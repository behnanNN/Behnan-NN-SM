## For GPU-node
salloc --time=02:00:00 -p gpu --gpus-per-node 1 --nodes 1 --ntasks 1 --cpus-per-task 16 --mem 128G

## For CPU-node
salloc --time=04:00:00 --nodes 1 --ntasks 1 --cpus-per-task 24 --mem 128G

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

cd45.obs['new_label_sample'] = cd45.obs['new_label'].astype('str') + "_"+ cd45.obs['sample'].astype('str')

## Removing lymphoid cells
cd45=cd45[~cd45.obs['new_label'].isin(['Lymphoid Lin.'])]

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

for i in tqdm(cd45.obs["new_label_sample"].unique().tolist()):
	bulk_ls.append(pseudobulk_sample(cd45,"new_label_sample",i,3,2))

bulk_dat = sc.concat(bulk_ls)
bulk_dat.to_df().to_csv("myeloid_pseudobulk.csv")

def mapper(val):
	return val[0:len(val)-2]

bulk_dat.obs['group'] = bulk_dat.obs.index.map(mapper)

metadata = pd.DataFrame(zip(bulk_dat.obs.index,bulk_dat.obs["group"])).rename(columns={0:"sample",1:"group"}).set_index("sample")
metadata.to_csv("myeloid_pseudobulk_metadata.csv")

## Setting up deseq2

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

bulk = pd.read_csv("myeloid_pseudobulk.csv", index_col=0)
metadata = pd.read_csv("myeloid_pseudobulk_metadata.csv",index_col=0)

dds = DeseqDataSet(counts=bulk,metadata=metadata,design_factors="group",refit_cooks=True,n_cpus=6)
sc.pp.filter_genes(dds,min_cells=1)

dds.deseq2()

dds.write('myeloid_pseudobulk_dds.h5ad')

def stats(object,grp1,grp2):
	stat_res = DeseqStats(object, contrast=('group',grp1,grp2))
	stat_res.summary()
	res = stat_res.results_df
	res = res[res.baseMean >= 5]
	#sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > np.log2(1.5))]
	#sigs = sigs.sort_values("log2FoldChange",ascending=False)
	return res

## Exporting Raw unfiltered lists to generate GSEA plots later on.

ctypes = dds.obs['group'].unique().tolist()
ctypes = [x.split('-')[0] for x in ctypes]
ctypes = list(set(ctypes))

master_analyzed_nocutoff = {}

for i,j in tqdm(zip(['gl_bm_vs_sham_bm','sb_bm_vs_sham_bm','gl_sm_vs_sham_sm','sb_sm_vs_sham_sm'],
	[('-gl-bm','-sham-bm'),('-sb-bm','-sham-bm'),('-gl-sk','-sham-sk'),('-sb-sk','-sham-sk')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_nocutoff[i]=temp_d

## Exporting no-cutoff files

for i in tqdm(list(master_analyzed_nocutoff)):
	os.makedirs('pydeseq2_analysis/'+i+'/nocutoff')
	path_to_save = 'pydeseq2_analysis/'+i+'/nocutoff/'
	temp_d = master_analyzed_nocutoff[i]
	for j in list(temp_d):
		temp_d[j].to_csv(path_to_save+j+'.csv')

## 1.2 Fold Change cutoff files

def stats_1_2(object,grp1,grp2):
	stat_res = DeseqStats(object, contrast=('group',grp1,grp2),lfc_null=np.log2(1.2))
	stat_res.summary()
	res = stat_res.results_df
	res = res[res.baseMean >= 5]
	sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > np.log2(1.2))]
	sigs = sigs.sort_values("log2FoldChange",ascending=False)
	return sigs

master_analyzed_1_2 = {}

for i,j in tqdm(zip(['gl_bm_vs_sham_bm','sb_bm_vs_sham_bm','gl_sm_vs_sham_sm','sb_sm_vs_sham_sm'],
	[('-gl-bm','-sham-bm'),('-sb-bm','-sham-bm'),('-gl-sk','-sham-sk'),('-sb-sk','-sham-sk')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats_1_2(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_1_2[i]=temp_d

## Exporting no-cutoff files

for i in tqdm(list(master_analyzed_1_2)):
	os.makedirs('pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff')
	path_to_save = 'pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff/'
	temp_d = master_analyzed_1_2[i]
	for j in list(temp_d):
		temp_d[j].to_csv(path_to_save+j+'.csv')


## Setting up plotting and analysis post DESEQ2

## Doing SM vs BM -> 3 Nov 2024

ctypes = dds.obs['group'].unique().tolist()
ctypes = [x.split('-')[0] for x in ctypes]
ctypes = list(set(ctypes))

master_analyzed_nocutoff = {}

for i,j in tqdm(zip(['sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm'],
	[('-sham-sk','-sham-bm'),('-gl-sk','-gl-bm'),('-sb-sk','-sb-bm')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_nocutoff[i]=temp_d

## Exporting no-cutoff files

for i in tqdm(list(master_analyzed_nocutoff)):
	os.makedirs('pydeseq2_analysis/'+i+'/nocutoff')
	path_to_save = 'pydeseq2_analysis/'+i+'/nocutoff/'
	temp_d = master_analyzed_nocutoff[i]
	for j in list(temp_d):
		temp_d[j].to_csv(path_to_save+j+'.csv')

## 1.2 Fold Change cutoff files

master_analyzed_1_2 = {}

for i,j in tqdm(zip(['sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm'],
	[('-sham-sk','-sham-bm'),('-gl-sk','-gl-bm'),('-sb-sk','-sb-bm')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats_1_2(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_1_2[i]=temp_d

## Exporting 1.2_fold_cutoff files

for i in tqdm(list(master_analyzed_1_2)):
	os.makedirs('pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff')
	path_to_save = 'pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff/'
	temp_d = master_analyzed_1_2[i]
	for j in list(temp_d):
		temp_d[j].to_csv(path_to_save+j+'.csv')

## Setting up plotting and analysis post DESEQ2