## For GPU-node
salloc --time=02:00:00 -p gpu --gpus-per-node 1 --nodes 1 --ntasks 1 --cpus-per-task 16 --mem 128G

## For CPU-node
salloc --time=04:00:00 --nodes 1 --ntasks 1 --cpus-per-task 24 --mem 128G

## Update on the previous analysis
## Midway, it was realized that sham bonemarrow that we had processed had something weird therefore needed to drop that
## Eventually, using TBM sample from one of the previously published studies


## Setting the directory

cd '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis'
conda activate scvi_tune

## Importing packages
import scanpy as sc
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import scvi
plt.rcParams['svg.fonttype'] = 'none'

sc.settings.figdir = './lymphoid_analysis/figures/'

## Importing cd45 file

cd45 = sc.read_h5ad('cd45.h5ad')

lym = cd45[cd45.obs['new_label'].isin(['Lymphoid Lin.'])].copy()
sc.pp.filter_genes(lym,min_cells=3)

lym.X = lym.layers['raw_counts'].copy()
sc.pp.normalize_total(lym, target_sum=1e+4)
sc.pp.log1p(lym)
lym.layers["log_counts"] = lym.X.copy()
sc.pp.highly_variable_genes(lym, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key='sample')
sc.pl.highly_variable_genes(lym,show=False,save='_hvg.png')
lym_hvg = lym[:,lym.var.highly_variable].copy()

## Basic filtering

sc.tl.pca(lym_hvg)
sc.pl.pca_variance_ratio(lym_hvg,n_pcs=50,show=False,log=True,save='_.png')
sc.pp.neighbors(lym_hvg, n_pcs=25)
sc.tl.umap(lym_hvg)
sc.tl.leiden(lym_hvg,flavor='igraph',resolution=0.7)

sc.pl.umap(lym_hvg,color=['sample'],show=False,size=15,save='_before_integration.png')

## umap splits
path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/lymphoid_analysis/figures/'

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(lym_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'before_integration_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')


## Setting up scvi to integrate

import torch
torch.set_float32_matmul_precision('high')

scvi.settings.num_threads = 8

scvi.model.SCVI.setup_anndata(lym_hvg, layer="raw_counts", batch_key="sample")
model = scvi.model.SCVI(lym_hvg, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

lym_hvg.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(lym_hvg, n_pcs=25,use_rep='X_scVI')
sc.tl.umap(lym_hvg)
sc.tl.leiden(lym_hvg,flavor='igraph',resolution=1)

## umap splits

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(lym_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'scvi_hvg_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

lym.obsm['X_umap'] = lym_hvg.obsm['X_umap']
lym.obs['leiden'] = lym_hvg.obs['leiden']

lym.write('lymphoid_analysis/lym.h5ad')
lym_hvg.write('lymphoid_analysis/lym_hvg.h5ad')


## Importing old annotations

dat_old = sc.read_h5ad('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/lym_hvg.h5ad')

old_d = dict(zip(dat_old.obs.index,dat_old.obs['lym_label']))
not_present = list(set(lym.obs.index) - set(dat_old.obs.index))

for i in not_present:
	old_d[i] = 'NA'

lym.obs['old_lym_label'] = lym.obs.index.map(old_d)
lym_hvg.obs['old_lym_label'] = lym_hvg.obs.index.map(old_d)

sc.pl.umap(lym_hvg,color=['leiden','old_lym_label'],legend_loc='on data',show=False,save='_old_label.png')

## Annotation

sc.pl.dotplot(lym,['Cd3g','Trac','Trdv4','Trgv2','Rora','Gata3','Cd4','Cd8a',
	'Klrb1c','Foxp3','Ncam1','Fcgr3','Rag1','Rag2','Dntt','Pax5','Chchd10',
	'Cd19','Ms4a1','H2-DMb2','Mzb1','Gypa','Cd79a','Kit','Spn','Cd34','Mki67',
	'Ly6g','Retnlg'],
	groupby='subcluster',show=False,save='_bcells.png',standard_scale='var')


sc.tl.leiden(lym_hvg,restrict_to=('leiden',['5']),resolution=0.2,key_added='subcluster')
lym_hvg.obs['leiden'] = lym_hvg.obs['subcluster']
sc.tl.leiden(lym_hvg,restrict_to=('leiden',['7']),resolution=0.2,key_added='subcluster')
lym_hvg.obs['leiden'] = lym_hvg.obs['subcluster']
lym.obs['leiden'] = lym_hvg.obs['leiden']

def annotator(val):
	if val=='9':
		return 'Antigen Presenting B-cells'
	if val in ['6','10']:
		return 'Mature B-cells'
	if val in ['3','11','13']:
		return 'Immature B-cells'
	if val=='7':
		return 'Pre B-cells'
	if val in ['2','1']:
		return 'Late Pro B-cells'
	if val=='12':
		return 'T-cells'
	if val=='15':
		return 'CD4 Tregs'
	if val=='8,0':
		return 'Early Pro B-cells'
	if val=='8,1':
		return 'Gamma Delta T-cells' 
	if val == '4,0':
		return 'NK-cells'
	if val in ['4,1','4,2']:
		return 'NKT-cells'
	if val == '14':
		return 'ILCs'
	else:
		return val

lym.obs['lym_label'] = lym.obs['leiden'].map(annotator)
lym_hvg.obs['lym_label'] = lym_hvg.obs['leiden'].map(annotator)

## splitting cls 8
sc.tl.leiden(lym_hvg,restrict_to=('leiden',['8']),resolution=0.15,key_added='subcluster')

sc.pl.umap(lym_hvg,color=['subcluster'],legend_loc='on data')

lym_hvg.obs['leiden'] = lym_hvg.obs['subcluster']
lym.obs['leiden'] = lym.obs['subcluster']

lym.obs['lym_label'] = lym.obs['leiden'].map(annotator)
lym_hvg.obs['lym_label'] = lym_hvg.obs['leiden'].map(annotator)

## splitting cls 4
sc.tl.leiden(lym_hvg,restrict_to=('leiden',['4']),resolution=0.15,key_added='subcluster')

## removing myeloid contaminant cells

lym = lym[~lym.obs['lym_label'].isin(['0','5'])]


## Trying to run paga with the new labels

sc.tl.paga(lym_hvg,'lym_label')
sc.pl.paga(lym_hvg,show=False)
sc.tl.umap(lym_hvg,init_pos='paga')

lym.obsm['X_umap'] = lym.obsm['X_umap']

lym.write('lym.h5ad')
lym_hvg.write('lym_hvg.h5ad')

------------------ Plotting and Figure export performed locally on ubuntu ---------------------
--------------------------- Script is on other file ------------------------------------------

## Relabeling clusters post removing myeloid cells.
## Old cluster information stored at 'lym_label_' in lym_hvg.obs

def annotator(val):
	if val=='9':
		return 'Antigen Presenting B-cells'
	if val in ['3','5','14']:
		return 'Mature B-cells'
	if val in ['8','6','7']:
		return 'Immature B-cells'
	if val=='11':
		return 'Pre B-cells'
	if val in ['10','4','0']:
		return 'Late Pro B-cells'
	if val=='2,1':
		return 'T-cells'
	if val=='15':
		return 'CD4 Tregs'
	if val in ['12','13,1']:
		return 'Early Pro B-cells'
	if val=='16':
		return 'Gamma Delta T-cells' 
	if val == '1':
		return 'NK-cells'
	if val in ['2,0','2,2']:
		return 'NKT-cells'
	if val == '13,0':
		return 'ILCs'
	else:
		return val