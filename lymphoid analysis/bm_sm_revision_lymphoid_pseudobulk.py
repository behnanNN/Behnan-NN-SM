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

lym=sc.read_h5ad('lymphoid_analysis/lym.h5ad')
lym_hvg=sc.read_h5ad('lymphoid_analysis/lym_hvg.h5ad')

lym.obs['lym_label_sample'] = lym.obs['lym_label'].astype('str') + "_"+ lym.obs['sample'].astype('str')

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

for i in tqdm(lym.obs["lym_label_sample"].unique().tolist()):
	bulk_ls.append(pseudobulk_sample(lym,"lym_label_sample",i,3,2))

bulk_dat = sc.concat(bulk_ls)
bulk_dat.to_df().to_csv("lymphoid_analysis/lymphoid_pseudobulk.csv")

def mapper(val):
	return val[0:len(val)-2]

bulk_dat.obs['group'] = bulk_dat.obs.index.map(mapper)

metadata = pd.DataFrame(zip(bulk_dat.obs.index,bulk_dat.obs["group"])).rename(columns={0:"sample",1:"group"}).set_index("sample")
metadata.to_csv("lymphoid_analysis/lymphoid_pseudobulk_metadata.csv")

## Setting up deseq2

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

bulk = pd.read_csv("lymphoid_analysis/lymphoid_pseudobulk.csv", index_col=0)
metadata = pd.read_csv("lymphoid_analysis/lymphoid_pseudobulk_metadata.csv",index_col=0)

dds = DeseqDataSet(counts=bulk,metadata=metadata,design_factors="group",refit_cooks=True,n_cpus=6)
sc.pp.filter_genes(dds,min_cells=1)

dds.deseq2()

dds.write('lymphoid_analysis/lymphoid_pseudobulk_dds.h5ad')

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

def ctype_mapper(val):
	if '-sham-' in val:
		return val.split('-sham-')[0]
	if '-gl-' in val:
		return val.split('-gl-')[0]
	if '-sb-' in val:
		return val.split('-sb-')[0]

ctypes = [ctype_mapper(x) for x in ctypes]
ctypes = list(set(ctypes))
ctypes.remove('Gamma Delta T-cells')
ctypes.remove('ILCs')

master_analyzed_nocutoff = {}

for i,j in tqdm(zip(['gl_bm_vs_sham_bm','sb_bm_vs_sham_bm','gl_sm_vs_sham_sm','sb_sm_vs_sham_sm','sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm'],
	[('-gl-bm','-sham-bm'),('-sb-bm','-sham-bm'),('-gl-sk','-sham-sk'),('-sb-sk','-sham-sk'),('-sham-sk','-sham-bm'),('-gl-sk','-gl-bm'),('-sb-sk','-sb-bm')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_nocutoff[i]=temp_d

## Exporting no-cutoff files

for i in tqdm(list(master_analyzed_nocutoff)):
	os.makedirs('lymphoid_analysis/pydeseq2_analysis/'+i+'/nocutoff')
	path_to_save = 'lymphoid_analysis/pydeseq2_analysis/'+i+'/nocutoff/'
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

for i,j in tqdm(zip(['gl_bm_vs_sham_bm','sb_bm_vs_sham_bm','gl_sm_vs_sham_sm','sb_sm_vs_sham_sm','sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm'],
	[('-gl-bm','-sham-bm'),('-sb-bm','-sham-bm'),('-gl-sk','-sham-sk'),('-sb-sk','-sham-sk'),('-sham-sk','-sham-bm'),('-gl-sk','-gl-bm'),('-sb-sk','-sb-bm')])):
	temp_d = {}
	for cellname in ctypes:
		sig = stats_1_2(dds,cellname+j[0],cellname+j[1])
		temp_d[cellname] = sig
	master_analyzed_1_2[i]=temp_d

## Exporting no-cutoff files

for i in tqdm(list(master_analyzed_1_2)):
	os.makedirs('lymphoid_analysis/pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff')
	path_to_save = 'lymphoid_analysis/pydeseq2_analysis/'+i+'/1.2_log2FC_0.05_pval_cutoff/'
	temp_d = master_analyzed_1_2[i]
	for j in list(temp_d):
		temp_d[j].to_csv(path_to_save+j+'.csv')


## Setting up plotting and analysis post DESEQ2