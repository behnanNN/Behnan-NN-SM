## For GPU-node
salloc --time=02:00:00 -p gpu --gpus-per-node 1 --nodes 1 --ntasks 1 --cpus-per-task 8 --mem 64G

## For CPU-node
salloc --time=04:00:00 --nodes 1 --ntasks 1 --cpus-per-task 24 --mem 128G

## Check the value of the following variables for the gpu-nodes. If CUDA_VISIBLE_DEVICES,
## has some predefined value, then please unset it else torch will not detect the gpu
## devices.

echo $CUDA_VISIBLE_DEVICES
echo $LD_LIBRARY_PATH

## Update on the previous analysis
## Sham BM in previous analysis was used from other published paper. Now we are replacing that sample in the new analysis

## Setting the directory

cd '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis'
conda activate scrnaseq

## Importing packages
import scvi
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import distinctipy
from tqdm import tqdm
plt.rcParams['svg.fonttype'] = 'none'
sc.set_figure_params(dpi_save=300, fontsize=10)

scvi.settings.num_threads = 8
sc.settings.figdir = 'preprocessing_figures/'
## Importing all the files 

## Using tdTomato mapped files

filepath = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/processed_fastqs/mapped_h5ads/tdTom_mapping/'

filenames = os.listdir(filepath)
filenames = [x for x in filenames if 'bm.h5ad' in x] + [x for x in filenames if 'sk.h5ad' in x]
filenames.remove('sham_bm.h5ad')

samples_d = {}

for i in tqdm(filenames):
	sample_name = i.split('.h5ad')[0]
	samples_d[sample_name] = sc.read_h5ad(filepath+i)

## Adding sample name as an observation to each sample

for i in tqdm(list(samples_d)):
	data=samples_d[i]
	data.obs['sample'] = i
	samples_d[i] = data

## Importing new sham bm sample
sham_bm = sc.read_10x_h5('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/flex_oct2024/BM2/outs/per_sample_outs/BM2/count/sample_filtered_feature_bc_matrix.h5')
sham_bm.obs['sample']='sham_bm'
sham_bm.var_names_make_unique()

samples_d['sham_bm'] = sham_bm

## Merging the samples into one sample

dat = sc.concat(list(samples_d.values()))
dat.obs_names_make_unique()
dat.var_names_make_unique()

## Basic filtering
sc.pp.filter_cells(dat,min_genes=200)
sc.pp.filter_genes(dat,min_cells=3)

## Mitochondrial QC
dat.var['mt'] = dat.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(dat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(dat,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,show=False,save='_QC.png')

>>> dat.shape
(47626, 16314)

dat = dat[dat.obs.pct_counts_mt < 10,:]

>>> dat.shape
(47185, 16314)

## Doublet Detection
sc.pp.scrublet(dat, batch_key="sample")

## Removing doublets

>>> dat.obs['predicted_doublet'].value_counts()
predicted_doublet
False    46548
True       637

dat = dat[dat.obs['predicted_doublet'].isin([False])]

## Normalization
dat.layers["raw_counts"] = dat.X.copy()
sc.pp.normalize_total(dat, target_sum=1e+4)
sc.pp.log1p(dat)
dat.layers["log_counts"] = dat.X.copy()
sc.pp.highly_variable_genes(dat, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key='sample')
sc.pl.highly_variable_genes(dat,show=False,save='_hvg.png')
dat_hvg = dat[:,dat.var.highly_variable].copy()

sc.tl.pca(dat_hvg)
sc.pl.pca_variance_ratio(dat_hvg,n_pcs=50,log=True,save='_.png')
sc.pp.neighbors(dat_hvg, n_pcs=30)
sc.tl.umap(dat_hvg)
sc.tl.leiden(dat_hvg,flavor='igraph')

sc.pl.umap(dat_hvg,color=['sample'],show=False,size=15,save='_before_integration.png')

## umap splits
path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/preprocessing_figures/'

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(dat_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'before_integration_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## Trying out harmony for integration

import scanpy.external as sce

sce.pp.harmony_integrate(dat_hvg,key='sample')
sc.pp.neighbors(dat_hvg, n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(dat_hvg)
sc.tl.leiden(dat_hvg,flavor='igraph')

sc.pl.umap(dat_hvg,color=['sample'],show=False,size=15,save='_after_integration_harmony.png')

## umap splits

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(dat_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'harmony_hvg_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## Setting up scvi to integrate

scvi.model.SCVI.setup_anndata(dat_hvg, layer="raw_counts", batch_key="sample")
model = scvi.model.SCVI(dat_hvg, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

dat_hvg.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(dat_hvg, n_pcs=30,use_rep='X_scVI')
sc.tl.umap(dat_hvg)
sc.tl.leiden(dat_hvg,flavor='igraph')

sc.pl.umap(dat_hvg,color=['sample'],show=False,size=15,save='_after_integration.png')

## umap splits

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(dat_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'scvi_hvg_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## Isolating only CD45+ cells

sc.pl.dotplot(dat,['Ptprc','Ly6g','Cd19','Cd3g','Kit','Cd34'],standard_scale='var',groupby='leiden',show=False,save='_cd45.png')

def cd45_annotator(val):
	if val in ['7','9','21','25','28']:
		return 'CD45neg'
	else:
		return 'CD45pos'

dat.obs['CD45'] = dat.obs['leiden'].map(cd45_annotator)
dat_hvg.obs['CD45'] = dat_hvg.obs['leiden'].map(cd45_annotator)

sc.pl.violin(dat,['Ptprc'],groupby='CD45',show=False,save='_cd45posneg.png')
sc.pl.dotplot(dat,['Ptprc'],groupby='CD45',show=False,save='_cd45posneg.png')

tmpls = []

for i in tqdm(dat.obs['sample'].unique().tolist()):
	d = dict(dat[dat.obs['sample'].isin([i])].obs['CD45'].value_counts())
	cd45_p = 100*d['CD45pos']/(d['CD45pos']+d['CD45neg'])
	cd45_n = 100*d['CD45neg']/(d['CD45pos']+d['CD45neg'])
	tmpls.append([d['CD45pos']+d['CD45neg'],d['CD45pos'],d['CD45neg'],cd45_p,cd45_n])

pd.DataFrame(tmpls,index=dat.obs['sample'].unique().tolist(),columns=['Total Cells','CD45 pos','CD45 neg','%CD45 pos','%CD45 neg']).to_csv('csvs/cd45_frequency.csv')

------------------------------------
## Interim saving of relevant files
dat.write('dat.h5ad')
dat_hvg.write('dat_hvg.h5ad')
------------------------------------

cd45 = dat[dat.obs['CD45'].isin(['CD45pos'])].copy()

cd45.X = cd45.layers['raw_counts'].copy()
sc.pp.normalize_total(cd45, target_sum=1e+4)
sc.pp.log1p(cd45)
cd45.layers["log_counts"] = cd45.X.copy()
sc.pp.highly_variable_genes(cd45, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key='sample')
sc.pl.highly_variable_genes(cd45,show=False,save='_hvg.png')
cd45_hvg = cd45[:,cd45.var.highly_variable].copy()

sc.tl.pca(cd45_hvg)
sc.pl.pca_variance_ratio(cd45_hvg,n_pcs=50,log=True,save='_.png')
sc.pp.neighbors(cd45_hvg, n_pcs=30)
sc.tl.umap(cd45_hvg)
sc.tl.leiden(cd45_hvg,flavor='igraph')
sc.pl.umap(cd45_hvg,color=['sample'],show=False,save='_before_integration.png')

## umap splits

path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/cd45_figures/'

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(cd45_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'before_integration_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')


## Moving to integration using scvi

scvi.model.SCVI.setup_anndata(cd45_hvg, layer="raw_counts", batch_key="sample")
model = scvi.model.SCVI(cd45_hvg, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

cd45_hvg.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(cd45_hvg, n_pcs=30,use_rep='X_scVI')
sc.tl.umap(cd45_hvg)
sc.tl.leiden(cd45_hvg,flavor='igraph',resolution=0.7)
sc.tl.paga(cd45_hvg,'leiden')
sc.pl.paga(cd45_hvg,show=False)
sc.tl.umap(cd45_hvg,init_pos='paga')

sc.pl.umap(cd45_hvg,color=['sample'],show=False,save='_after_integration.png')

cd45.obs['leiden'] = cd45_hvg.obs['leiden']
cd45.obsm['X_umap'] = cd45_hvg.obsm['X_umap']

## umap splits

fig,ax=plt.subplots(figsize=(12,6),nrows=2,ncols=3)

for i,j in zip(['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm'],ax.ravel()):
	sc.pl.umap(cd45_hvg,color=['sample'],groups=[i],size=15,show=False,ax=j)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'after_integration_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')


## Mapping old clusters to the new space

dat_old = sc.read_h5ad('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/cd45_hvg.h5ad')

old_d = dict(zip(dat_old.obs.index,dat_old.obs['mye_label']))

not_present = list(set(cd45.obs.index) - set(dat_old.obs.index))

for i in not_present:
	old_d[i] = 'NA'

cd45.obs['old_label'] = cd45.obs.index.map(old_d)
cd45_hvg.obs['old_label'] = cd45_hvg.obs.index.map(old_d)

sc.pl.umap(cd45,color=['leiden','old_label'],legend_loc='on data',show=False,save='_old_label.png')

## umap split of old_label

old_label = ['Erythro. Prog.', 'Monocytes', 'Pre. Neutro.', 'MEPs', 'GMPs',
 'Prol. Monocytes', 'STHSCs', 'Lymphoid Lin.', 'DCs', 'Prol. Macrophages',
  'Transit. Neutro. 1', 'Mature Neutro.', 'Hemato. Prog.', 'Macrophages',
  'Transit. Neutro. 2', 'Acp5+ Macrophages']

fig,ax=plt.subplots(figsize=(16,16),nrows=4,ncols=4)

for i,j in zip(old_label,ax.ravel()):
	sc.pl.umap(cd45_hvg,color=['old_label'],groups=[i],size=15,show=False,ax=j,legend_loc='none')
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'old_label_UMAP_split.png',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

def cls_annot(val):
	if val in ['10','9']:
		return 'Mature Neutro.'
	if val == '4':
		return 'Transit. Neutro. 1'
	if val in ['5','11,0']:
		return 'Transit. Neutro. 2'
	if val in ['12','13','17','18','19','2,3']:
		return 'Lymphoid Lin.'
	if val == '0':
		return 'Erythro. Prog.'
	if val == '7':
		return 'STHSCs'
	if val == '15':
		return 'Hemato. Prog.'
	if val in ['14','11,2','1']:
		return 'Macrophages'
	if val == '16':
		return 'Acp5+ Macrophages'
	if val in ['6,0','6,1']:
		return 'Prol. Macrophages'
	if val == '6,2':
		return 'Monocytes'	
	if val == '8':
		return 'DCs'
	if val in ['3','2,0']:
		return 'Pre. Neutro.'
	if val in ['2,1','11,1']:
		return 'GMPs'
	if val == '2,2':
		return 'MEPs'
	else:
		return val

cd45.obs['new_label'] = cd45.obs['leiden'].map(cls_annot)
cd45_hvg.obs['new_label'] = cd45_hvg.obs['leiden'].map(cls_annot)

## splitting cls 2
sc.tl.leiden(cd45_hvg,restrict_to=('leiden',['2']),resolution=0.15,key_added='subcluster')

sc.pl.umap(cd45_hvg,color=['subcluster'],legend_loc='on data')

cd45_hvg.obs['leiden'] = cd45_hvg.obs['subcluster']
cd45.obs['leiden'] = cd45.obs['subcluster']

cd45.obs['new_label'] = cd45.obs['subcluster'].map(cls_annot)
cd45_hvg.obs['new_label'] = cd45_hvg.obs['subcluster'].map(cls_annot)

## splitting cls 11
sc.tl.leiden(cd45_hvg,restrict_to=('leiden',['11']),resolution=0.15,key_added='subcluster')

## splitting cls 6
sc.tl.leiden(cd45_hvg,restrict_to=('leiden',['6']),resolution=0.15,key_added='subcluster')

## Trying to run paga with the new labels

sc.tl.paga(cd45_hvg,'new_label')
sc.pl.paga(cd45_hvg,show=False)
sc.tl.umap(cd45_hvg,init_pos='paga')

cd45.obsm['X_umap'] = cd45_hvg.obsm['X_umap']

cd45.write('cd45.h5ad')
cd45_hvg.write('cd45_hvg.h5ad')

------------------ Plotting and Figure export performed locally on ubuntu ---------------------
--------------------------- Script is on other file ------------------------------------------