import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import distinctipy
from tqdm import tqdm
from matplotlib import patches
plt.rcParams['svg.fonttype'] = 'none'
sc.set_figure_params(dpi_save=300, fontsize=10)

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

plot_order = ['Lymphoid Lin.','Neutrophils','Monocytes','Macrophages','DCs','Hemato./Erythro Prog.','HSCs']

mye_df = pd.read_csv('facs_toplot/myeloid_facs_toplot.csv',index_col=0)
mye_sk_df = mye_df[['sham_sk','gl_sk','sb_sk']].loc[plot_order].copy()
mye_bm_df = mye_df[['sham_bm','gl_bm','sb_bm']].loc[plot_order].copy()

## Myeloid colors_d

import pickle

with open('myeloid_color_d.pickle','rb') as handle:
	colors_d=pickle.load(handle)

mye_facs_colors_d = {}

mye_facs_colors_d['Lymphoid Lin.'] = colors_d['Lymphoid Lin.']
mye_facs_colors_d['Neutrophils'] = colors_d['Mature Neutro.']
mye_facs_colors_d['Monocytes'] = colors_d['Monocytes']
mye_facs_colors_d['Macrophages'] = colors_d['Macrophages']
mye_facs_colors_d['DCs'] = colors_d['DCs']
mye_facs_colors_d['Hemato./Erythro Prog.'] = colors_d['Erythro. Prog.']
mye_facs_colors_d['HSCs'] = colors_d['STHSCs']

## Plotting Sankey Diagram

axx=pysankey(mye_sk_df.T,ncols=mye_sk_df.T.shape[0],color_d=mye_facs_colors_d)
plt.savefig("facs_toplot/_sk_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/_sk_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

axx=pysankey(mye_bm_df.T,ncols=mye_bm_df.T.shape[0],color_d=mye_facs_colors_d)
plt.savefig("facs_toplot/_bm_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/_bm_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")


## Similar Plotting for lymphoid

lym_plot_order = ['B-cells','CD4 T-cells','CD4 Tregs','CD8 T-cells','Gamma Delta T-cells',
'DN- T-cells','NK/NKT Cells']

lym_df = pd.read_csv('facs_toplot/lymphoid_facs_toplot.csv',index_col=0)
lym_sk_df = lym_df[['sham_sk','gl_sk','sb_sk']].loc[lym_plot_order].copy()
lym_bm_df = lym_df[['sham_bm','gl_bm','sb_bm']].loc[lym_plot_order].copy()


lym_hvg=sc.read_h5ad('lymphoid_analysis/lym_hvg.h5ad')
lym_colors_d=dict(zip(lym_hvg.obs['lym_label'].unique(),lym_hvg.uns['lym_label_colors']))

lym_facs_colors_d = {}

lym_facs_colors_d['B-cells'] = lym_colors_d['Immature B-cells']
lym_facs_colors_d['CD4 T-cells'] = 'grey'
lym_facs_colors_d['CD4 Tregs'] = lym_colors_d['CD4 Tregs']
lym_facs_colors_d['CD8 T-cells'] = lym_colors_d['CD8 T-cells']
lym_facs_colors_d['Gamma Delta T-cells'] = lym_colors_d['Gamma Delta T-cells']
lym_facs_colors_d['DN- T-cells'] = 'maroon'
lym_facs_colors_d['NK/NKT Cells'] = lym_colors_d['NKT-cells']

fig,ax=plt.subplots()
lym_sk_df.T.plot.bar(stacked=True,color=lym_facs_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_sk_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_sk_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

fig,ax=plt.subplots()
lym_bm_df.T.plot.bar(stacked=True,color=lym_facs_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_bm_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_bm_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

## Plotting lymphoid out of cd45 for scrnaseq and then facs

plot_order = ['Early Pro B-cells','Late Pro B-cells','Pre B-cells','Immature B-cells',
'Mature B-cells','Antigen Presenting B-cells','Plasma B-cells','CD8 T-cells','CD4 Tregs',
'NK-cells','NKT-cells','Gamma Delta T-cells','Myeloid']

lym_cd45_scrna_df = pd.read_csv('facs_toplot/lymphoid_scrnaseq_toplot_cd45.csv',index_col=0)

lym_cd45sk_scrna_df = lym_cd45_scrna_df[['sham_sk','gl_sk','sb_sk']]
lym_cd45bm_scrna_df = lym_cd45_scrna_df[['sham_bm','gl_bm','sb_bm']]

lym_hvg=sc.read_h5ad('lymphoid_analysis/lym_hvg.h5ad')
lym_colors_d=dict(zip(lym_hvg.obs['lym_label'].unique(),lym_hvg.uns['lym_label_colors']))

lym_colors_d['Myeloid'] = 'orange'

fig,ax=plt.subplots()
lym_cd45sk_scrna_df.T.plot.bar(stacked=True,color=lym_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_CD45_scRNAseq_sk_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_CD45_scRNAseq_sk_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

fig,ax=plt.subplots()
lym_cd45bm_scrna_df.T.plot.bar(stacked=True,color=lym_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_CD45_scRNAseq_bm_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_CD45_scRNAseq_bm_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")


## Doing same cd45 plotting for flow data

lym_facs_colors_d['Myeloid'] = 'orange'

lym_plot_order = ['B-cells','CD4 T-cells','CD4 Tregs','CD8 T-cells','Gamma Delta T-cells',
'DN- T-cells','NK/NKT Cells','Myeloid']

lym_cd45_facs_df = pd.read_csv('facs_toplot/lymphoid_facs_toplot_cd45.csv',index_col=0)

lym_cd45sk_facs_df = lym_cd45_facs_df[['sham_sk','gl_sk','sb_sk']]
lym_cd45bm_facs_df = lym_cd45_facs_df[['sham_bm','gl_bm','sb_bm']]

fig,ax=plt.subplots()
lym_cd45sk_facs_df.T.plot.bar(stacked=True,color=lym_facs_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_CD45_facs_sk_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_CD45_facs_sk_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")

fig,ax=plt.subplots()
lym_cd45bm_facs_df.T.plot.bar(stacked=True,color=lym_facs_colors_d,yticks=[0,10,20,30,40,50,60,70,80,90,100]).legend(
	loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.savefig("facs_toplot/lym_CD45_facs_bm_df"+".png",dpi=300, bbox_inches="tight")
plt.savefig("facs_toplot/lym_CD45_facs_bm_df"+".svg",dpi=300, bbox_inches="tight")
plt.clf()
plt.close("all")