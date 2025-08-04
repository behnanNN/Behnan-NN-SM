import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import distinctipy
from matplotlib import patches
plt.rcParams['svg.fonttype'] = 'none'
sc.set_figure_params(dpi_save=300, fontsize=10)
sc.settings.figdir = './lymphoid_analysis/figures/'

lym = sc.read_h5ad('lymphoid_analysis/lym.h5ad')
lym_hvg = sc.read_h5ad('lymphoid_analysis/lym_hvg.h5ad')

## Defining colors dict and saving it too

colors = distinctipy.get_colors(12)
colors_d = dict(zip(lym_hvg.obs['lym_label'].unique(),colors))

import pickle

with open('lymphoid_analysis/lymphoid_color_d.pickle', 'wb') as handle:
    pickle.dump(colors_d, handle, protocol=pickle.HIGHEST_PROTOCOL)

## Exporting cluster annotation UMAP
sc.pl.umap(lym_hvg,color=['lym_label'],palette=colors_d,show=False,size=8,save='_cluster_annotation.png')
sc.pl.umap(lym_hvg,color=['lym_label'],palette=colors_d,show=False,size=8,save='_cluster_annotation.svg')

## Now making split umaps

path_to_save='lymphoid_analysis/figures/'

fig,ax = plt.subplots(figsize=(9,6),nrows=2,ncols=3)

grp_names = ['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm']

for i,j in zip(grp_names,ax.ravel()):
	pol = sc.pl.umap(lym_hvg,show=False,ax=j,frameon=False,size=15,legend_loc=None)
	sc.pl.umap(lym_hvg[lym_hvg.obs['sample'].isin([i])],color=['lym_label'],palette=colors_d,size=15,frameon=False,show=False,ax=j,legend_loc=None)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'_umap_splits.png',bbox_inches='tight',dpi=300)
plt.savefig(path_to_save+'_umap_splits.svg',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## Frequency Dataframe

plot_order = ['Early Pro B-cells','Late Pro B-cells','Pre B-cells','Immature B-cells',
'Mature B-cells','Antigen Presenting B-cells','T-cells','CD4 Tregs',
'NK-cells','NKT-cells','Gamma Delta T-cells','ILCs']

out_ls = []
grp_ls = ['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm']
cls_ls = lym_hvg.obs["lym_label"].unique().tolist()
total_count_d = dict(lym_hvg.obs["sample"].value_counts())

for i in tqdm(grp_ls):
	subs_dat = lym_hvg[lym_hvg.obs["sample"].isin([i])]
	count_d = dict(subs_dat.obs["lym_label"].value_counts())
	tmp_ls=[]
	for j in tqdm(cls_ls):
		if j not in list(count_d):
			tmp_ls.append(0)
		else:
			tmp_ls.append(100*count_d[j]/total_count_d[i])
	out_ls.append(tmp_ls)

cls_df=pd.DataFrame(out_ls,index=grp_ls,columns=cls_ls).T
cls_df=cls_df.loc[plot_order]
cls_df.to_csv("lymphoid_analysis/csvs/lymphoid_cluster_frequencies.csv")

sk_df = cls_df[['sham_sk','gl_sk','sb_sk']]
sk_df=sk_df.loc[~(sk_df==0).all(axis=1)]
bm_df = cls_df[['sham_bm','gl_bm','sb_bm']]
bm_df=bm_df.loc[~(sk_df==0).all(axis=1)]

## Making counts dataframe as well. 
out_ls_counts = []

for i in tqdm(grp_ls):
	subs_dat = lym_hvg[lym_hvg.obs["sample"].isin([i])]
	count_d = dict(subs_dat.obs["lym_label"].value_counts())
	tmp_ls=[]
	for j in tqdm(cls_ls):
		if j not in list(count_d):
			tmp_ls.append(0)
		else:
			tmp_ls.append(count_d[j])
	out_ls_counts.append(tmp_ls)

counts_df=pd.DataFrame(out_ls_counts,index=grp_ls,columns=cls_ls).T
counts_df=counts_df.loc[plot_order]
counts_df.to_csv("lymphoid_analysis/csvs/lymphoid_cluster_counts.csv")

## Making Sankey diagrams

def pysankey(dff,ncols,color_d,figsize=(5,6)):
	fig,ax = plt.subplots(figsize=figsize)
	pol=dff.plot.bar(stacked=True,ax=ax,color=color_d)
	patches_ = ax.patches
	same_patches_ = [patches_[x:x+ncols] for x in range(0,len(patches_),ncols)]
	for i,col in zip(same_patches_,[color_d[x] for x in dff.columns]):
		for j in range(0,len(i)-1):
			first_patch = i[j]
			next_patch = i[j+1]
			p1 = first_patch.get_corners()[1]
			p2 = next_patch.get_corners()[0]		
			p3 = next_patch.get_corners()[3]
			p4 = first_patch.get_corners()[2]			
			ax.add_patch(patches.Polygon([p1,p2,p3,p4],linewidth=0,facecolor=col,alpha=0.7))
	return ax

axx=pysankey(sk_df.T,ncols=sk_df.T.shape[0],color_d=colors_d)
plt.savefig('lymphoid_analysis/figures/sankey_trial_sk.png',dpi=300)
plt.savefig('lymphoid_analysis/figures/sankey_trial_sk.svg',dpi=300)
plt.clf()
plt.close('all')

axx=pysankey(bm_df.T,ncols=sk_df.T.shape[0],color_d=colors_d)
plt.savefig('lymphoid_analysis/figures/sankey_trial_bm.png',dpi=300)
plt.savefig('lymphoid_analysis/figures/sankey_trial_bm.svg',dpi=300)
plt.clf()
plt.close('all')


@@ START LATER

## Making an annotation dotplot (TO BE DONE LATER)

cat_order = ['Early Pro B-cells','Late Pro B-cells','Pre B-cells','Immature B-cells',
'Mature B-cells','Antigen Presenting B-cells','T-cells','CD4 Tregs',
'NK-cells','NKT-cells','Gamma Delta T-cells','ILCs']


lym.obs['lym_label'] = lym.obs['lym_label'].cat.reorder_categories(cat_order)
lym_hvg.obs['lym_label'] = lym_hvg.obs['lym_label'].cat.reorder_categories(cat_order)

markers_dd = ['Rag1', 'Rag2', 'Dntt', 'Chchd10', 'Il7r', 'Mki67', 'Cd19', 'Ms4a1', 'H2-Aa', 'H2-Eb1', 'Cd3g', 'Cd8a', 'Cd4', 'Trac', 'Ctla4', 'Foxp3', 'Klrb1c', 'Nkg7', 'Cd3d', 'Rora','Gata3']


## For the purpose of fig, need to change the dimensions of the plotting
## Compressed version for plotting

fig,ax = plt.subplots(figsize=(8,3.5))
sc.pl.dotplot(lym,markers_dd,groupby='lym_label',standard_scale='var',cmap='PRGn_r',show=False,ax=ax)
plt.tight_layout()
plt.savefig(path_to_save+'lymphoid_final_annotation_compressed.png',bbox_inches='tight',dpi=300)
plt.savefig(path_to_save+'lymphoid_final_annotation_compressed.svg',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## END ##