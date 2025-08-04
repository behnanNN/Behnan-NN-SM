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

cd45 = sc.read_h5ad('cd45.h5ad')
cd45_hvg = sc.read_h5ad('cd45_hvg.h5ad')

## Defining colors dict and saving it too

colors = distinctipy.get_colors(15)
colors_d = dict(zip(cd45_hvg.obs['new_label'].unique(),colors))

import pickle

with open('myeloid_color_d.pickle', 'wb') as handle:
    pickle.dump(colors_d, handle, protocol=pickle.HIGHEST_PROTOCOL)

## Exporting cluster annotation UMAP
sc.pl.umap(cd45_hvg,color=['new_label'],palette=colors_d,show=False,size=8,save='_cluster_annotation.png')
sc.pl.umap(cd45_hvg,color=['new_label'],palette=colors_d,show=False,size=8,save='_cluster_annotation.svg')

## Now making split umaps

path_to_save='cd45_figures/'

fig,ax = plt.subplots(figsize=(9,6),nrows=2,ncols=3)

grp_names = ['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm']

for i,j in zip(grp_names,ax.ravel()):
	pol = sc.pl.umap(cd45_hvg,show=False,ax=j,frameon=False,size=8,legend_loc=None)
	sc.pl.umap(cd45_hvg[cd45_hvg.obs['sample'].isin([i])],color=['new_label'],palette=colors_d,size=8,frameon=False,show=False,ax=j,legend_loc=None)
	j.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'_umap_splits.png',bbox_inches='tight',dpi=300)
plt.savefig(path_to_save+'_umap_splits.svg',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## Frequency Dataframe

plot_order = ['Lymphoid Lin.','Mature Neutro.','Transit. Neutro. 1','Transit. Neutro. 2','Pre. Neutro.',
'Monocytes','Macrophages','Acp5+ Macrophages','Prol. Macrophages','DCs','Erythro. Prog.',
'Hemato. Prog.','GMPs','MEPs','STHSCs']

out_ls = []
grp_ls = ['sham_sk','gl_sk','sb_sk','sham_bm','gl_bm','sb_bm']
cls_ls = cd45_hvg.obs["new_label"].unique().tolist()
total_count_d = dict(cd45_hvg.obs["sample"].value_counts())

for i in tqdm(grp_ls):
	subs_dat = cd45_hvg[cd45_hvg.obs["sample"].isin([i])]
	count_d = dict(subs_dat.obs["new_label"].value_counts())
	tmp_ls=[]
	for j in tqdm(cls_ls):
		if j not in list(count_d):
			tmp_ls.append(0)
		else:
			tmp_ls.append(100*count_d[j]/total_count_d[i])
	out_ls.append(tmp_ls)

cls_df=pd.DataFrame(out_ls,index=grp_ls,columns=cls_ls).T
cls_df=cls_df.loc[plot_order]
cls_df.to_csv("csvs/myeloid_cluster_frequencies.csv")

sk_df = cls_df[['sham_sk','gl_sk','sb_sk']]
sk_df=sk_df.loc[~(sk_df==0).all(axis=1)]
bm_df = cls_df[['sham_bm','gl_bm','sb_bm']]
bm_df=bm_df.loc[~(sk_df==0).all(axis=1)]

## Making counts dataframe as well. 
out_ls_counts = []

for i in tqdm(grp_ls):
	subs_dat = cd45_hvg[cd45_hvg.obs["sample"].isin([i])]
	count_d = dict(subs_dat.obs["new_label"].value_counts())
	tmp_ls=[]
	for j in tqdm(cls_ls):
		if j not in list(count_d):
			tmp_ls.append(0)
		else:
			tmp_ls.append(count_d[j])
	out_ls_counts.append(tmp_ls)

counts_df=pd.DataFrame(out_ls_counts,index=grp_ls,columns=cls_ls).T
counts_df=counts_df.loc[plot_order]
counts_df.to_csv("csvs/myeloid_cluster_counts.csv")

## Now plotting cluster frequencies stacked bar chart
## (Will make alluvial plot later on google colab R session)

fig,ax=plt.subplots()
sk_df.T.plot.bar(stacked=True,color=colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig(path_to_save+"_sk_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig(path_to_save+"_sk_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

fig,ax=plt.subplots()
bm_df.T.plot.bar(stacked=True,color=colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig(path_to_save+"_bm_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig(path_to_save+"_bm_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

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
plt.savefig('cd45_figures/sankey_trial_sk.png',dpi=300)
plt.savefig('cd45_figures/sankey_trial_sk.svg',dpi=300)
plt.clf()
plt.close('all')

axx=pysankey(bm_df.T,ncols=bm_df.T.shape[0],color_d=colors_d)
plt.savefig('cd45_figures/sankey_trial_bm.png',dpi=300)
plt.savefig('cd45_figures/sankey_trial_bm.svg',dpi=300)
plt.clf()
plt.close('all')

## Making an annotation dotplot

with open('myeloid_color_d.pickle','rb') as handle:
	colors_d=pickle.load(handle)


cat_order = ['STHSCs','GMPs','MEPs','Hemato. Prog.','Erythro. Prog.','Monocytes','DCs',
'Macrophages','Prol. Macrophages','Acp5+ Macrophages','Pre. Neutro.','Transit. Neutro. 1',
'Transit. Neutro. 2','Mature Neutro.','Lymphoid Lin.']

cd45.obs['new_label'] = cd45.obs['new_label'].cat.reorder_categories(cat_order)
cd45_hvg.obs['new_label'] = cd45_hvg.obs['new_label'].cat.reorder_categories(cat_order)

markers_d = {'':['Ptprc'],'Stem/Stem Prog.':['Kit','Cd34','Flt3','Cd27','Ly6a','Fcgr3'],
'Erythro Prog.':['Hba-a1','Hba-a2','Gypa'],'Monocytes':['Ly6c2','Ccr2','Mki67'],
'DCs':['Itgax','Ccr9','Clec9a'],
'Macrophages':['Adgre1','Cd68','Mki67','Mrc1','Pf4','Acp5','H2-Aa','H2-Ab1'],
'Neutrophils':['Chil3','Camp','Ngp','Ly6g','Retnlg','Cxcr2'],'Lymphoid':['Cd3e','Cd19','Klrb1c']}

markers_dd = ['Ptprc','Kit','Cd34','Flt3','Cd27','Ly6a','Fcgr3','Hba-a1','Hba-a2','Gypa','Ccr2',
'Mki67','Itgax','Ccr9','Clec9a','Cd68','Mki67','Mrc1','Pf4','Acp5','H2-Aa','H2-Ab1','Camp',
'Ngp','Ly6g','Retnlg','Cxcr2','Cd3e','Cd19','Klrb1c']

sc.pl.dotplot(cd45,markers_d,groupby='mye_label',standard_scale='var',cmap='Spectral_r',show=False,save='_myeloid_final_annotation_Acp5.png')
sc.pl.dotplot(cd45,markers_d,groupby='mye_label',standard_scale='var',cmap='Spectral_r',show=False,save='_myeloid_final_annotation.svg')

## For the purpose of fig, need to change the dimensions of the plotting
## Compressed version for plotting

fig,ax = plt.subplots(figsize=(8,6))
sc.pl.dotplot(cd45,markers_dd,groupby='new_label',standard_scale='var',cmap='RdBu_r',show=False,ax=ax)
plt.tight_layout()
plt.savefig(path_to_save+'_myeloid_final_annotation_compressed.png',bbox_inches='tight',dpi=300)
plt.savefig(path_to_save+'_myeloid_final_annotation_compressed.svg',bbox_inches='tight',dpi=300)
plt.clf()
plt.close('all')

## END ##