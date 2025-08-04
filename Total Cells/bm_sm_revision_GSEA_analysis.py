import gseapy as gp
from gseapy import dotplot
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from tqdm import tqdm
plt.rcParams['svg.fonttype'] = 'none'

## mouse to human mapper

from gseapy import Biomart
bm = Biomart()
# note the dataset and attribute names are different
m2h = bm.query(dataset='mmusculus_gene_ensembl',
               attributes=['ensembl_gene_id','external_gene_name',
                           'hsapiens_homolog_ensembl_gene',
                           'hsapiens_homolog_associated_gene_name'])

## Saving this for future conversions
m2h.to_csv('mouse_to_human_conversion.csv')

m2h_d = dict(zip(m2h['external_gene_name'],m2h['hsapiens_homolog_associated_gene_name']))

## Doing GSEA analysis for all the comparisons
dir_names = os.listdir('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/pydeseq2_analysis')
# dir_names = [x for x in dir_names if '.png' not in x]
# dir_names = [x for x in dir_names if '.svg' not in x]
# dir_names.remove('summary_plots_overlapping_genes')

file_list=os.listdir('/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/pydeseq2_analysis/sb_sm_vs_sham_sm/nocutoff')
genesets = 'GO_Biological_Process_2023'

common_path='/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/pydeseq2_analysis/'

for i in tqdm(dir_names):
	os.makedirs(common_path+i+'/GSEA_analysis/GO_Biological_Process/plots')
	path_to_save = common_path+i+'/GSEA_analysis/GO_Biological_Process/'
	for j in tqdm(file_list):
		df=pd.read_csv(common_path+i+'/nocutoff/'+j,index_col=0)
		df['h_gene'] = df.index.map(m2h_d)
		df=df.reset_index().set_index('h_gene')
		pre_res = gp.prerank(rnk=df["log2FoldChange"],gene_sets='GO_Biological_Process_2023',threads=20,min_size=5,
                     max_size=1000,permutation_num=1000,outdir=None,seed=6,verbose=True)
		title = j.split(".csv")[0]
		dff = pre_res.res2d
		dff = dff[dff['FDR q-val']<=0.25]
		if dff.shape[0]!=0:
			dff.to_csv(path_to_save+title+' '+i+"_GO Biological Process.csv")	
			fig,ax= plt.subplots()
			ax = dotplot(pre_res.res2d,column="FDR q-val",title=title+"_"+i,cmap=plt.cm.viridis,
				size=3,figsize=(15,10), cutoff=0.25, show_ring=False)
			plt.tight_layout()
			plt.savefig(path_to_save+'plots/'+title+"_"+i+"_GO Biological Process.png")
			plt.clf()
			plt.close("all")
		else:
			pass
			dff.to_csv(path_to_save+title+' '+i+"_GO Biological Process.csv")
	
	
## For each cell-type and across conditions, retrieve total differential pathways as well as directions	

pathways_d = {}
pathways_summary_ = []

for i in tqdm(dir_names):
	path_to_look = common_path+i+'/GSEA_analysis/GO_Biological_Process/'
	tmpd = {}
	tmpls = []
	cell_types = os.listdir(path_to_look)
	cell_types.remove('plots')
	for j in cell_types:
		title=j.split(' '+i)[0]
		dff=pd.read_csv(path_to_look+j,index_col=0)
		tmpd[title] = dff
		if dff.shape[0]!=0:
			dff_up = dff[dff['NES']>0]
			dff_dn = dff[dff['NES']<0]
			tmpls.append([dff.shape[0],dff_up.shape[0],dff_dn.shape[0]])
		else:
			tmpls.append([0,0,0])
	tmpdf = pd.DataFrame(tmpls,index=[x.split(''+i)[0] for x in cell_types],columns=[i+'_Total',i+'_Up',i+'_Dn'])
	pathways_summary_.append(tmpdf)	
	pathways_d[i]=tmpd

pathways_total_df = pd.concat(pathways_summary_,axis=1)

pathways_all = pathways_total_df[[x+'_Total' for x in dir_names]]
pathways_up = pathways_total_df[[x+'_Up' for x in dir_names]]
pathways_dn = pathways_total_df[[x+'_Dn' for x in dir_names]]

pathways_all.to_csv('csvs/GSEA_analysis/summary_total_differential_pathways_across_groups.csv')
pathways_up.to_csv('csvs/GSEA_analysis/summary_Upregulated_pathways_across_groups.csv')
pathways_dn.to_csv('csvs/GSEA_analysis/summary_Downregulated_pathways_across_groups.csv')

## Let's plot first part of the figure before moving forward. 

cat_order = ['STHSCs','GMPs','MEPs','Hemato. Prog.','Erythro. Prog.','Monocytes',
'DCs','Macrophages','Prol. Macrophages',
'Acp5+ Macrophages','Pre. Neutro.','Transit. Neutro. 1','Transit. Neutro. 2','Mature Neutro.']

cat_order.reverse()

for i in tqdm(dir_names):
	tmpdf=pathways_total_df[[i+x for x in ['_Up','_Dn']]].copy()
	tmpdf[i+'_Dn'] = tmpdf[i+'_Dn'].multiply(-1)
	tmpdf.index = tmpdf.index.map(lambda x:x[0:len(x)-1])
	tmpdf=tmpdf.loc[cat_order]
	tmpdf.to_csv('csvs/GSEA_analysis/'+'Up_down_pathways_'+i+'.csv')
	fig,ax = plt.subplots(figsize=(5,3))
	tmpdf[i+'_Up'].plot.barh(ax=ax,color='Red',width=0.9,xlim=(-550,600),xticks=list(range(-550,600,100)))
	tmpdf[i+'_Dn'].plot.barh(ax=ax,color='Blue',width=0.9,xlim=(-550,600),xticks=list(range(-550,600,100)))
	for container in ax.containers:
		ax.bar_label(container,padding=5)
	ax.vlines(x=0,ymin=-1,ymax=15,color='black')
	ax.title.set_text(i)
	plt.xticks(rotation=90)
	plt.tight_layout()
	plt.savefig('gsea_plots/'+i+'_up_down_plot.png',dpi=300)
	plt.savefig('gsea_plots/'+i+'_up_down_plot.svg',dpi=300)
	plt.clf()
	plt.close('all')

## Now trying to pull distinct and common pathways for sm analysis first

combined_d = {}

for i in tqdm(list(pathways_d)):
	tmpd = pathways_d[i]
	all_pathways = [x for k in list(tmpd) for x in tmpd[k]['Term'].unique().tolist()]
	all_pathways = list(set(all_pathways))
	tmpls = []
	for cell_types in list(tmpd):
		tmpdf = tmpd[cell_types]
		nes_d = dict(zip(tmpd[cell_types]['Term'],tmpd[cell_types]['NES']))
		nes_d.update(dict.fromkeys(list(set(all_pathways)-set(list(nes_d))),0))
		tmpls.append([nes_d[x] for x in all_pathways])
	outdf=pd.DataFrame(tmpls,index=list(tmpd),columns=all_pathways)
	combined_d[i] = outdf

## Using sc.tl.rank_genes_groups_df to get top markers for each group

adatas_d = {}

for i in tqdm(list(combined_d)):
	dff=combined_d[i]
	adata=anndata.AnnData(pd.concat([dff,dff,dff],axis=0))
	adata.obs_names_make_unique()
	adata.obs['Group'] = adata.obs.index.map(lambda x:x.split('-')[0])
	adata.obs['Group'] = adata.obs['Group'].astype('category')
	sc.tl.rank_genes_groups(adata,'Group',method='t-test_overestim_var')
	adatas_d[i] = adata

for i in tqdm(list(adatas_d)):
	adata=adatas_d[i]
	fig,ax=plt.subplots(figsize=(20,10))
	sc.pl.rank_genes_groups_dotplot(adata,show=False,dendrogram=False,ax=ax)
	plt.tight_layout()
	plt.savefig('gsea_plots/pathways_'+i+'.png',bbox_inches='tight',dpi=300)
	plt.clf()
	plt.close('all')

## Saving for tomorrow

for i in tqdm(list(adatas_d)):
	adatas_d[i].write('pathways_adata/pathways_'+i+'.h5ad')


## Starting 10 June 2024 Analysis

## Make a barplot showing pathways commonly up or down in response to both tumors. Make this
## for both SM and BM

glsm_d = pathways_d['gl_sm_vs_sham_sm']
sbsm_d = pathways_d['sb_sm_vs_sham_sm']

sm_shared_up_d = {}
sm_shared_dn_d = {}
sm_shared_summary_ls = []

for i in tqdm(cat_order):
	gldf = glsm_d[i].copy()
	sbdf = sbsm_d[i].copy()
	if gldf.shape[0]!=0 and sbdf.shape[0]!=0:
		gldf_up = gldf[gldf['NES']>0].set_index('Term')
		gldf_dn = gldf[gldf['NES']<0].set_index('Term')
		sbdf_up = sbdf[sbdf['NES']>0].set_index('Term')
		sbdf_dn = sbdf[sbdf['NES']<0].set_index('Term')
		shared_up = list(set(gldf_up.index).intersection(set(sbdf_up.index)))
		shared_dn = list(set(gldf_dn.index).intersection(set(sbdf_dn.index))) 
		shared_up_df = pd.concat([gldf_up.loc[shared_up][['NES']].rename(columns={'NES':'NES_GL'}),
			sbdf_up.loc[shared_up][['NES']].rename(columns={'NES':'NES_SB'})],axis=1)
		shared_up_df=shared_up_df.sort_values('NES_GL',ascending=False)
		shared_dn_df = pd.concat([gldf_dn.loc[shared_dn][['NES']].rename(columns={'NES':'NES_GL'}),
			sbdf_dn.loc[shared_dn][['NES']].rename(columns={'NES':'NES_SB'})],axis=1)
		shared_dn_df=shared_dn_df.sort_values('NES_GL')
		sm_shared_up_d[i] = shared_up_df
		sm_shared_dn_d[i] = shared_dn_df
		sm_shared_summary_ls.append([len(shared_up),-1*len(shared_dn)])
	else:
		sm_shared_summary_ls.append([0,0])

sm_shared_summary_df = pd.DataFrame(sm_shared_summary_ls,index=cat_order,columns=['Shared_Upregulated_SM','Shared_Downregulated_SM'])
sm_shared_summary_df.to_csv('csvs/GSEA_shared_pathways/summary_shared_sm.csv')

## Plotting similar bar plot

fig,ax = plt.subplots(figsize=(4,3))
sm_shared_summary_df['Shared_Upregulated_SM'].plot.barh(ax=ax,color='Red',width=0.9,xlim=(-100,150),xticks=list(range(-100,200,50)))
sm_shared_summary_df['Shared_Downregulated_SM'].plot.barh(ax=ax,color='Blue',width=0.9,xlim=(-100,150),xticks=list(range(-100,200,50)))
for container in ax.containers:
	ax.bar_label(container,padding=5)

ax.vlines(x=0,ymin=-1,ymax=15,color='black')
ax.title.set_text('Tumor vs ShamSM Shared Pathways')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('gsea_plots/shared_up_down/'+'SM_shared_up_down_pathways.png',dpi=300)
plt.savefig('gsea_plots/shared_up_down/'+'SM_shared_up_down_pathways.svg',dpi=300)
plt.clf()
plt.close('all')


path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/csvs/GSEA_shared_pathways/sm/'

for i in tqdm(cat_order):
	if i in list(sm_shared_up_d) and list(sm_shared_dn_d):
		dfup=sm_shared_up_d[i]
		dfdn=sm_shared_dn_d[i]
		dfup.to_csv(path_to_save+'upregulated/'+i+'_Shared_Upregulated_SM.csv')
		dfdn.to_csv(path_to_save+'downregulated/'+i+'_Shared_Downregulated_SM.csv')
	else:
		pass

## Combining shared up and shared dn across cell types

for i in tqdm(list(sm_shared_up_d)):
	dfup=sm_shared_up_d[i]
	dfdn=sm_shared_dn_d[i]
	dfup=dfup.reset_index()
	dfdn=dfdn.reset_index()
	dfup['Cell Type'] = i
	dfdn['Cell Type'] = i
	sm_shared_up_d[i] = dfup.set_index('Cell Type')
	sm_shared_dn_d[i] = dfdn.set_index('Cell Type')

combined_sm_shared_up_df=pd.concat(list(sm_shared_up_d.values()),axis=0)
combined_sm_shared_dn_df=pd.concat(list(sm_shared_dn_d.values()),axis=0)

combined_sm_shared_up_df.to_csv('csvs/GSEA_shared_pathways/sm_shared_upregulated_pathways_all_celltypes.csv')
combined_sm_shared_dn_df.to_csv('csvs/GSEA_shared_pathways/sm_shared_downregulated_pathways_all_celltypes.csv')


## Repeating similar analysis for bone marrow

glbm_d = pathways_d['gl_bm_vs_sham_bm']
sbbm_d = pathways_d['sb_bm_vs_sham_bm']

bm_shared_up_d = {}
bm_shared_dn_d = {}
bm_shared_summary_ls = []

for i in tqdm(cat_order):
	gldf = glbm_d[i].copy()
	sbdf = sbbm_d[i].copy()
	if gldf.shape[0]!=0 and sbdf.shape[0]!=0:
		gldf_up = gldf[gldf['NES']>0].set_index('Term')
		gldf_dn = gldf[gldf['NES']<0].set_index('Term')
		sbdf_up = sbdf[sbdf['NES']>0].set_index('Term')
		sbdf_dn = sbdf[sbdf['NES']<0].set_index('Term')
		shared_up = list(set(gldf_up.index).intersection(set(sbdf_up.index)))
		shared_dn = list(set(gldf_dn.index).intersection(set(sbdf_dn.index))) 
		shared_up_df = pd.concat([gldf_up.loc[shared_up][['NES']].rename(columns={'NES':'NES_GL'}),
			sbdf_up.loc[shared_up][['NES']].rename(columns={'NES':'NES_SB'})],axis=1)
		shared_up_df=shared_up_df.sort_values('NES_GL',ascending=False)
		shared_dn_df = pd.concat([gldf_dn.loc[shared_dn][['NES']].rename(columns={'NES':'NES_GL'}),
			sbdf_dn.loc[shared_dn][['NES']].rename(columns={'NES':'NES_SB'})],axis=1)
		shared_dn_df=shared_dn_df.sort_values('NES_GL')
		bm_shared_up_d[i] = shared_up_df
		bm_shared_dn_d[i] = shared_dn_df
		bm_shared_summary_ls.append([len(shared_up),-1*len(shared_dn)])
	else:
		bm_shared_summary_ls.append([0,0])

bm_shared_summary_df = pd.DataFrame(bm_shared_summary_ls,index=cat_order,columns=['Shared_Upregulated_BM','Shared_Downregulated_BM'])
bm_shared_summary_df.to_csv('csvs/GSEA_shared_pathways/summary_shared_bm.csv')

## Plotting similar bar plot

fig,ax = plt.subplots(figsize=(4,3))
bm_shared_summary_df['Shared_Upregulated_BM'].plot.barh(ax=ax,color='Red',width=0.9,xlim=(-300,50),xticks=list(range(-300,100,50)))
bm_shared_summary_df['Shared_Downregulated_BM'].plot.barh(ax=ax,color='Blue',width=0.9,xlim=(-300,50),xticks=list(range(-300,100,50)))
for container in ax.containers:
	ax.bar_label(container,padding=5)

ax.vlines(x=0,ymin=-1,ymax=15,color='black')
ax.title.set_text('Tumor vs ShamBM Shared Pathways')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('gsea_plots/shared_up_down/'+'BM_shared_up_down_pathways.png',dpi=300)
plt.savefig('gsea_plots/shared_up_down/'+'BM_shared_up_down_pathways.svg',dpi=300)
plt.clf()
plt.close('all')


path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/csvs/GSEA_shared_pathways/bm/'

for i in tqdm(cat_order):
	if i in list(bm_shared_up_d) and list(bm_shared_dn_d):
		dfup=bm_shared_up_d[i]
		dfdn=bm_shared_dn_d[i]
		dfup.to_csv(path_to_save+'upregulated/'+i+'_Shared_Upregulated_BM.csv')
		dfdn.to_csv(path_to_save+'downregulated/'+i+'_Shared_Downregulated_BM.csv')
	else:
		pass

## Combining shared up and shared dn across cell types

for i in tqdm(list(bm_shared_up_d)):
	dfup=bm_shared_up_d[i]
	dfdn=bm_shared_dn_d[i]
	dfup=dfup.reset_index()
	dfdn=dfdn.reset_index()
	dfup['Cell Type'] = i
	dfdn['Cell Type'] = i
	bm_shared_up_d[i] = dfup.set_index('Cell Type')
	bm_shared_dn_d[i] = dfdn.set_index('Cell Type')

combined_bm_shared_up_df=pd.concat(list(bm_shared_up_d.values()),axis=0)
combined_bm_shared_dn_df=pd.concat(list(bm_shared_dn_d.values()),axis=0)

combined_bm_shared_up_df.to_csv('csvs/GSEA_shared_pathways/bm_shared_upregulated_pathways_all_celltypes.csv')
combined_bm_shared_dn_df.to_csv('csvs/GSEA_shared_pathways/bm_shared_downregulated_pathways_all_celltypes.csv')


## Looking to combine most common pathways affected across cell types

most_up_sbsm = []
most_dn_sbsm = []
most_up_glsm = []
most_dn_glsm = []

for i in tqdm(list(combined_d['sb_sm_vs_sham_sm'])):
	dff=combined_d['sb_sm_vs_sham_sm'][i]
	dffup=dff[dff>0]
	dffdn=dff[dff<0]
	most_up_sbsm.append([i,dffup.shape[0]])
	most_dn_sbsm.append([i,dffdn.shape[0]])

for i in tqdm(list(combined_d['gl_sm_vs_sham_sm'])):
	dff=combined_d['gl_sm_vs_sham_sm'][i]
	dffup=dff[dff>0]
	dffdn=dff[dff<0]
	most_up_glsm.append([i,dffup.shape[0]])
	most_dn_glsm.append([i,dffdn.shape[0]])

sbsm_mostupdf = pd.DataFrame(most_up_sbsm).rename(columns={0:'Term',1:'N'}).set_index('Term')
sbsm_mostdndf = pd.DataFrame(most_dn_sbsm).rename(columns={0:'Term',1:'N'}).set_index('Term')
glsm_mostupdf = pd.DataFrame(most_up_glsm).rename(columns={0:'Term',1:'N'}).set_index('Term')
glsm_mostdndf = pd.DataFrame(most_dn_glsm).rename(columns={0:'Term',1:'N'}).set_index('Term')

sbsm_mostupdf.sort_values('N',ascending=False,inplace=True)
sbsm_mostdndf.sort_values('N',inplace=True)
glsm_mostupdf.sort_values('N',ascending=False,inplace=True)
glsm_mostdndf.sort_values('N',inplace=True)

shared_top_up = set(sbsm_mostupdf.iloc[0:50].index).intersection(glsm_mostupdf.iloc[0:50].index)
shared_top_dn = set(sbsm_mostdndf.iloc[0:50].index).intersection(glsm_mostdndf.iloc[0:50].index)


---------------------
@@ Start Working Here
---------------------

## Performing semantic clustering first in R and then using that to define the list of pathways to be plotted

smup_df = combined_sm_shared_up_df.reset_index()['Term']
smdn_df = combined_sm_shared_dn_df.reset_index()['Term']
bmup_df = combined_bm_shared_up_df.reset_index()['Term']
bmdn_df = combined_bm_shared_dn_df.reset_index()['Term']

smup_df = pd.DataFrame(smup_df.map(lambda x:'GO'+x.split('GO')[1].split(')')[0])).rename(columns={'Term':'sm_up'})
smdn_df = pd.DataFrame(smdn_df.map(lambda x:'GO'+x.split('GO')[1].split(')')[0])).rename(columns={'Term':'sm_dn'})
bmup_df = pd.DataFrame(bmup_df.map(lambda x:'GO'+x.split('GO')[1].split(')')[0])).rename(columns={'Term':'bm_up'})
bmdn_df = pd.DataFrame(bmdn_df.map(lambda x:'GO'+x.split('GO')[1].split(')')[0])).rename(columns={'Term':'bm_dn'})

smup_df.to_csv('csvs/For_semantic_clustering/sm_up_forclustering.csv')
smdn_df.to_csv('csvs/For_semantic_clustering/sm_dn_forclustering.csv')
bmup_df.to_csv('csvs/For_semantic_clustering/bm_up_forclustering.csv')
bmdn_df.to_csv('csvs/For_semantic_clustering/bm_dn_forclustering.csv')

---------
@Moving to R for semantic clustering
---------

## Importing semantic clustering information from R csv files and plotting 

## go_id mapper

smup_df = pd.read_csv('SM_upregulated_tolook.csv',index_col=0)
smup_df['GO_ids'] = smup_df['0'].map(lambda x:'GO'+x.split('GO')[1].split(')')[0])
smup_semantic_d = dict(zip(smup_df['GO_ids'],smup_df['0']))

smup_semantic_df = pd.read_csv('GO_semantic_clustering/semantic_cluster_SM_UP.csv',index_col=0)
smup_semantic_df['Term'] = smup_semantic_df['id'].map(smup_semantic_d)
smup_semantic_df.groupby('cluster').agg({'Term':lambda x: list(x)})

dff_toplot = pd.concat([combined_d['sb_sm_vs_sham_sm'][smup_semantic_df.iloc[0]['Term']],combined_d['gl_sm_vs_sham_sm'][smup_semantic_df.iloc[0]['Term']]],axis=0).reset_index().groupby('index').mean()

tmpls = []

for i in dff_toplot.index:
	tmpls.append(dff_toplot.loc[i].mean().item())

dff_toplot_mean=pd.DataFrame(tmpls,index=dff_toplot.index).T


## Useless


## Plotting select pathways for SM and BM

sm_shared = pd.read_csv('SM_shared_toplot.csv')
bm_shared = pd.read_csv('BM_shared_toplot.csv')

## defining sm_exlcude to remove pathways that are no more in the list

sm_exclude = ['Regulation Of Interleukin-10 Production (GO:0032653)', 'Platelet Aggregation (GO:0070527)', 
'Negative Regulation Of Macrophage Activation (GO:0043031)', 'Negative Regulation Of Interleukin-1 Production (GO:0032692)', 
'Mast Cell Activation Involved In Immune Response (GO:0002279)', 'Macrophage Activation (GO:0042116)', 
'T Cell Receptor Signaling Pathway (GO:0050852)', 'Negative Regulation Of Angiogenesis (GO:0016525)', 
'Negative Regulation Of Cytokine Production (GO:0001818)', 'Regulation Of Adipose Tissue Development (GO:1904177)', 
'Epidermal Growth Factor Receptor Signaling Pathway (GO:0007173)', 'Response To Glucocorticoid (GO:0051384)', 
'Regulation Of Blood Coagulation (GO:0030193)', 'Positive Regulation Of Type II Interferon Production (GO:0032729)', 
'Dendritic Cell Chemotaxis (GO:0002407)', 'Positive Regulation Of Axonogenesis (GO:0050772)']


sm_shared = sm_shared[~sm_shared['Term'].isin(sm_exclude)]

dff_toplot_sm = pd.concat([combined_d['sb_sm_vs_sham_sm'][sm_shared['Term']],combined_d['gl_sm_vs_sham_sm'][sm_shared['Term']]],axis=0).reset_index().groupby('index').mean().T
dff_toplot_sm['Terms'] = dff_toplot_sm.index.map(lambda x:x.split(' (GO')[0])
dff_toplot_sm = dff_toplot_sm.set_index('Terms')
dff_toplot_sm['for_sort']=dff_toplot_sm.sum(axis=1)
dff_toplot_sm=dff_toplot_sm.sort_values('for_sort')
del dff_toplot_sm['for_sort']

sns.set(rc = {'figure.figsize':(9.6,3)},font_scale=0.6)
sns.heatmap(dff_toplot_sm.T.loc[cat_order],cmap='RdBu_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('gsea_plots/SM_select_pathways.svg')
plt.clf()
plt.close('all')


## for bm

bm_shared = bm_shared[~bm_shared['Term'].isin(bm_exclude)]

dff_toplot_bm = pd.concat([combined_d['sb_bm_vs_sham_bm'][bm_shared['Term']],combined_d['gl_bm_vs_sham_bm'][bm_shared['Term']]],axis=0).reset_index().groupby('index').mean().T
dff_toplot_bm['Terms'] = dff_toplot_bm.index.map(lambda x:x.split(' (GO')[0])
dff_toplot_bm = dff_toplot_bm.set_index('Terms')
dff_toplot_bm['for_sort']=dff_toplot_bm.sum(axis=1)
dff_toplot_bm=dff_toplot_bm.sort_values('for_sort')
del dff_toplot_bm['for_sort']

sns.set(rc = {'figure.figsize':(8,3)},font_scale=0.6)
sns.heatmap(dff_toplot_bm.T.loc[cat_order],cmap='RdBu_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('gsea_plots/BM_select_pathways.svg')
plt.clf()
plt.close('all')


## Overlap of bone marrow and skull marrow shared pathways

smup_df = pd.read_csv('SM_upregulated_tolook.csv',index_col=0)
smdn_df = pd.read_csv('SM_downregulated_tolook.csv',index_col=0)
bmup_df = pd.read_csv('BM_upregulated_tolook.csv',index_col=0)
bmdn_df = pd.read_csv('BM_downregulated_tolook.csv',index_col=0)

dict_for_venn = {'SM_UP':set(smup_df['0']),'SM_DN':set(smdn_df['0']),'BM_UP':set(bmup_df['0']),'BM_DN':set(bmdn_df['0'])}

## Now looking for pathways going up and down in SM

opp_pathways = list(dict_for_venn['SM_UP'].intersection(dict_for_venn['BM_DN'])) + list(dict_for_venn['SM_DN'].intersection(dict_for_venn['BM_UP']))
opp_pathways_filter_df=pd.read_csv('opp_pathways_selected.csv')
## Removing big name pathways
opp_pathways.remove('Positive Regulation Of Adaptive Immune Response Based On Somatic Recombination Of Immune Receptors Built From Immunoglobulin Superfamily Domains (GO:0002824)')
opp_pathways.remove('Maturation Of SSU-rRNA From Tricistronic rRNA Transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA) (GO:0000462)')

dff_opp_sm = pd.concat([combined_d['sb_sm_vs_sham_sm'][opp_pathways_filter_df['Terms']],combined_d['gl_sm_vs_sham_sm'][opp_pathways_filter_df['Terms']]],axis=0).reset_index().groupby('index').mean().T
dff_opp_bm = pd.concat([combined_d['sb_bm_vs_sham_bm'][opp_pathways_filter_df['Terms']],combined_d['gl_bm_vs_sham_bm'][opp_pathways_filter_df['Terms']]],axis=0).reset_index().groupby('index').mean().T

dff_opp_sm.columns=[x+'_SM' for x in dff_opp_sm.columns]
dff_opp_bm.columns=[x+'_BM' for x in dff_opp_bm.columns]
dff_opp_bm['for_sort']=dff_opp_bm.sum(axis=1)
dff_opp_bm=dff_opp_bm.sort_values('for_sort',ascending=False)
del dff_opp_bm['for_sort']
dff_opp_sm=dff_opp_sm.loc[dff_opp_bm.index]

dff_opp = pd.concat([dff_opp_sm,dff_opp_bm],axis=1)
dff_opp['Terms'] = dff_opp.index.map(lambda x:x.split(' (GO')[0])
dff_opp = dff_opp.set_index('Terms')

opp_cat_order = []
for i in cat_order:
	opp_cat_order.append(i+'_SM')
	opp_cat_order.append(i+'_BM')

sns.set(rc = {'figure.figsize':(5.6,6)},font_scale=0.6)
pol=sns.heatmap(dff_opp.T.loc[opp_cat_order],cmap='RdBu_r',xticklabels=True)
#plt.tight_layout()
plt.savefig('gsea_plots/opposite_pathways.svg')
plt.clf()
plt.close('all')



## Loading GSEA output files

pathways_d = {}

for i in os.listdir( '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis/pydeseq2_analysis/sb_sm_vs_sham_sm/go_biological_process/'):
	title = i.split(".csv")[0]
	pathways_d[title]=pd.read_csv(path_to_save+i,index_col=0)

out = []

for i in list(pathways_d):
	df = pathways_d[i]
	df_up = df[df['NES']>0]
	df_dn = df[df['NES']<0]
	out.append([i,df_up.shape[0],df_dn.shape[0]])

pd.DataFrame(out).set_index(0).rename(columns={1:"Upregulated_Pathways",2:"Downregulated_Pathways"}).to_csv("Pathways_summary.csv")

term_ls = []

for i in list(pathways_d):
	df = pathways_d[i]
	term_ls = term_ls + df["Term"].to_list()

term_ls = list(set(term_ls))

## 370 unique terms !!

## Total upregultaed and downregulated pathways

up_ls = []
dn_ls = []

for i in list(pathways_d):
	df = pathways_d[i]
	dfup = df[df['NES']>0]
	dfdn = df[df['NES']<0]
	up_ls = up_ls + dfup["Term"].to_list()
	dn_ls = dn_ls + dfdn["Term"].to_list()

up_ls = list(set(up_ls))
dn_ls = list(set(dn_ls))

out_df = []

for i in tqdm(term_ls):
	res = []
	for j in list(pathways_d):
		if i in pathways_d[j]["Term"].to_list():
			if pathways_d[j][pathways_d[j]['Term'].isin([i])]['NES'].item()>0:
				res.append('Up')
			else:
				res.append('Dn')
		else:
			res.append("No")
	out_df.append(res)

pathways_df = pd.DataFrame(out_df,index=term_ls,columns=list(pathways_d))

## Seeing which are the top ranked pathways among most cell types

out_dfdf = []

for i in tqdm(pathways_df.index.tolist()):
	dff = pathways_df.loc[i]
	tmpd = dict(pathways_df.loc[i].value_counts())
	not_present = list(set(['No','Up','Dn'])-set(list(tmpd)))
	if len(not_present)>0:
		for k in not_present:
			tmpd[k] = 0
	else:
		pass
	out_dfdf.append([tmpd['Up'],tmpd['Dn'],tmpd['No']])
	
pathways_summary_df = pd.DataFrame(out_dfdf,index=pathways_df.index.tolist(),columns=['Up','Dn','No'])
pathways_summary_df = pathways_summary_df.sort_values("Up",ascending=False)
pathways_summary_df_c = pathways_summary_df.sort_values("Dn",ascending=False)

## Plotting top 25 pathways this way

out_ls = []

for j in tqdm(pathways_summary_df.index.tolist()):
	tmpls = []
	for i in tqdm(list(pathways_d)):
		dff=pathways_d[i].set_index('Term')
		if j in dff.index.tolist():
			tmpls.append(dff.loc[j]['NES'].item())
		else:
			tmpls.append(0)
	out_ls.append(tmpls)

comb_df = pd.DataFrame(out_ls,index=pathways_summary_df.index.tolist(),columns=list(pathways_d))

fig,ax = plt.subplots(figsize=(8,10))
sns.heatmap(comb_df,cmap='bwr',ax=ax)
plt.tight_layout()
plt.savefig('trial_pathways_.png')

## Plotting barplot for 1st part of the figure as dicsussed

df['Downregulated_Pathways'] = -1*df['Downregulated_Pathways']

fig,ax = plt.subplots(figsize=(5,4))

df['Upregulated_Pathways'].plot.barh(ax=ax,color='Red',width=0.8,xlim=(-400,600))
df['Downregulated_Pathways'].plot.barh(ax=ax,color='Blue',width=0.8,xlim=(-400,600))

for container in ax.containers:
	ax.bar_label(container,padding=5)

ax.vlines(x=0,ymin=-1,ymax=10,color='black')
plt.tight_layout()
plt.savefig('temp_up_down_path.svg')


## Venn diagram for Tumor SM vs Sham SM

from venn import venn

dict_for_venn = {'GL261':set(combined_d['gl_sm_vs_sham_sm'].columns),
'SB28':set(combined_d['sb_sm_vs_sham_sm'].columns)}

venn(dict_for_venn,cmap=['red','green'])
plt.tight_layout()
plt.savefig('gsea_plots/venn_tumor_sm_vs_sham_sm.png',dpi=300)
plt.savefig('gsea_plots/venn_tumor_sm_vs_sham_sm.svg',dpi=300)
plt.clf()
plt.close('all')

## Venn diagram for Tumor BM vs Sham BM

from venn import venn

dict_for_venn = {'GL261':set(combined_d['gl_bm_vs_sham_bm'].columns),
'SB28':set(combined_d['sb_bm_vs_sham_bm'].columns)}

venn(dict_for_venn,cmap=['red','green'])
plt.tight_layout()
plt.savefig('gsea_plots/venn_tumor_bm_vs_sham_bm.png',dpi=300)
plt.savefig('gsea_plots/venn_tumor_bm_vs_sham_bm.svg',dpi=300)
plt.clf()
plt.close('all')