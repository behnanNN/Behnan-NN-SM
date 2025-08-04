import gseapy as gp
from gseapy import dotplot
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import seaborn as sns
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
dir_names = os.listdir('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis')
dir_names.remove('summary_plots_overlapping_genes')

file_list=os.listdir('/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/sb_sm_vs_sham_sm/nocutoff')
genesets = 'GO_Biological_Process_2023'

common_path='/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/'

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

pathways_all.to_csv('lymphoid_analysis/csvs/GSEA_analysis/summary_total_differential_pathways_across_groups.csv')
pathways_up.to_csv('lymphoid_analysis/csvs/GSEA_analysis/summary_Upregulated_pathways_across_groups.csv')
pathways_dn.to_csv('lymphoid_analysis/csvs/GSEA_analysis/summary_Downregulated_pathways_across_groups.csv')

## Let's plot first part of the figure before moving forward. 

cat_order = ['Early Pro B-cells','Late Pro B-cells','Pre B-cells','Immature B-cells',
'Mature B-cells','Antigen Presenting B-cells','CD8 T-cells','CD4 Tregs',
'NK-cells','NKT-cells']

cat_order.reverse()

for i in tqdm(dir_names):
	tmpdf=pathways_total_df[[i+x for x in ['_Up','_Dn']]].copy()
	tmpdf[i+'_Dn'] = tmpdf[i+'_Dn'].multiply(-1)
	tmpdf.index = tmpdf.index.map(lambda x:x[0:len(x)-1])
	tmpdf=tmpdf.loc[cat_order]
	tmpdf.to_csv('lymphoid_analysis/csvs/GSEA_analysis/'+'Up_down_pathways_'+i+'.csv')
	fig,ax = plt.subplots(figsize=(5,3))
	tmpdf[i+'_Up'].plot.barh(ax=ax,color='Red',width=0.9,xlim=(-550,600),xticks=list(range(-550,600,100)))
	tmpdf[i+'_Dn'].plot.barh(ax=ax,color='Blue',width=0.9,xlim=(-550,600),xticks=list(range(-550,600,100)))
	for container in ax.containers:
		ax.bar_label(container,padding=5)
	ax.vlines(x=0,ymin=-1,ymax=15,color='black')
	ax.title.set_text(i)
	plt.xticks(rotation=90)
	plt.tight_layout()
	plt.savefig('lymphoid_analysis/gsea_plots/'+i+'_up_down_plot.png',dpi=300)
	plt.savefig('lymphoid_analysis/gsea_plots/'+i+'_up_down_plot.svg',dpi=300)
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
	plt.savefig('pathways_'+i+'.png',bbox_inches='tight',dpi=300)
	plt.clf()
	plt.close('all')

## Saving for tomorrow

for i in tqdm(list(adatas_d)):
	adatas_d[i].write('lymphoid_analysis/pathways_adata/pathways_'+i+'.h5ad')


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
sm_shared_summary_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/summary_shared_sm.csv')

## Plotting similar bar plot

fig,ax = plt.subplots(figsize=(4,3))
sm_shared_summary_df['Shared_Upregulated_SM'].plot.barh(ax=ax,color='Purple',width=0.9,xlim=(-75,50),xticks=list(range(-75,75,25)))
sm_shared_summary_df['Shared_Downregulated_SM'].plot.barh(ax=ax,color='Green',width=0.9,xlim=(-75,50),xticks=list(range(-75,75,25)))
for container in ax.containers:
	ax.bar_label(container,padding=5)

ax.vlines(x=0,ymin=-1,ymax=15,color='black')
ax.title.set_text('Tumor vs ShamSM Shared Pathways')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/shared_up_down/'+'SM_shared_up_down_pathways.png',dpi=300)
plt.savefig('lymphoid_analysis/gsea_plots/shared_up_down/'+'SM_shared_up_down_pathways.svg',dpi=300)
plt.clf()
plt.close('all')


path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/csvs/GSEA_shared_pathways/sm/'

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

combined_sm_shared_up_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sm_shared_upregulated_pathways_all_celltypes.csv')
combined_sm_shared_dn_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sm_shared_downregulated_pathways_all_celltypes.csv')


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
bm_shared_summary_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/summary_shared_bm.csv')

## Plotting similar bar plot

fig,ax = plt.subplots(figsize=(4,3))
bm_shared_summary_df['Shared_Upregulated_BM'].plot.barh(ax=ax,color='Purple',width=0.9,xlim=(-75,50),xticks=list(range(-75,75,25)))
bm_shared_summary_df['Shared_Downregulated_BM'].plot.barh(ax=ax,color='Green',width=0.9,xlim=(-75,50),xticks=list(range(-75,75,25)))
for container in ax.containers:
	ax.bar_label(container,padding=5)

ax.vlines(x=0,ymin=-1,ymax=15,color='black')
ax.title.set_text('Tumor vs ShamBM Shared Pathways')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/shared_up_down/'+'BM_shared_up_down_pathways.png',dpi=300)
plt.savefig('lymphoid_analysis/gsea_plots/shared_up_down/'+'BM_shared_up_down_pathways.svg',dpi=300)
plt.clf()
plt.close('all')


path_to_save = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/csvs/GSEA_shared_pathways/bm/'

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

combined_bm_shared_up_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/bm_shared_upregulated_pathways_all_celltypes.csv')
combined_bm_shared_dn_df.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/bm_shared_downregulated_pathways_all_celltypes.csv')


## June 14 Analysis

## importing semantic clustering

gldf=pd.DataFrame(combined_d['gl_sm_vs_gl_bm'].columns)
sbdf=pd.DataFrame(combined_d['sb_sm_vs_sb_bm'].columns)

gldf['go_ids'] = gldf[0].map(lambda x:'GO'+x.split('GO')[1].split(')')[0])
sbdf['go_ids'] = sbdf[0].map(lambda x:'GO'+x.split('GO')[1].split(')')[0])

gldf.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/gl_sm_vs_bm_combined_pathways.csv')
sbdf.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sb_sm_vs_bm_combined_pathways.csv')


gl_semantic_df = pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/semantic_cluster_GL_SM_vs_BM.csv',index_col=0)
sb_semantic_df = pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/semantic_cluster_SB_SM_vs_BM.csv',index_col=0)

gl_semantic_d = dict(zip(gl_semantic_df['id'],gl_semantic_df['cluster']))
sb_semantic_d = dict(zip(sb_semantic_df['id'],sb_semantic_df['cluster']))

gldf['cluster'] = gldf['go_ids'].map(gl_semantic_d)
gldf=gldf.sort_values('cluster')
sbdf['cluster'] = sbdf['go_ids'].map(sb_semantic_d)
sbdf=sbdf.sort_values('cluster')

gldf.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/gl_sm_vs_bm_combined_pathways.csv')
sbdf.to_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sb_sm_vs_bm_combined_pathways.csv')


## Now plotting selected pathways

sm_shared = pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sm_shared_toplot.csv')
bm_shared = pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/bm_shared_toplot.csv')

dff_toplot_sm = pd.concat([combined_d['sb_sm_vs_sham_sm'][sm_shared['sm']],combined_d['gl_sm_vs_sham_sm'][sm_shared['sm']]],axis=0).reset_index().groupby('index').mean().T
dff_toplot_sm['Terms'] = dff_toplot_sm.index.map(lambda x:x.split(' (GO')[0])
dff_toplot_sm = dff_toplot_sm.set_index('Terms')
dff_toplot_sm['for_sort']=dff_toplot_sm.sum(axis=1)
dff_toplot_sm=dff_toplot_sm.sort_values('for_sort')
del dff_toplot_sm['for_sort']

sns.set(rc = {'figure.figsize':(4.2,2)},font_scale=0.6)
sns.heatmap(dff_toplot_sm.T.loc[cat_order],cmap='PRGn_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/SM_select_pathways.svg')
plt.clf()
plt.close('all')

## Alternate plotting vertically (IGNORE OTHERWISE)

sns.set(rc = {'figure.figsize':(2,4.2)},font_scale=0.6)
sns.heatmap(dff_toplot_sm[cat_order],cmap='PRGn_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/SM_select_pathways_alternate.svg')
plt.clf()
plt.close('all')

## for bm

dff_toplot_bm = pd.concat([combined_d['sb_bm_vs_sham_bm'][bm_shared['bm']],combined_d['gl_bm_vs_sham_bm'][bm_shared['bm']]],axis=0).reset_index().groupby('index').mean().T
dff_toplot_bm['Terms'] = dff_toplot_bm.index.map(lambda x:x.split(' (GO')[0])
dff_toplot_bm = dff_toplot_bm.set_index('Terms')
dff_toplot_bm['for_sort']=dff_toplot_bm.sum(axis=1)
dff_toplot_bm=dff_toplot_bm.sort_values('for_sort')
del dff_toplot_bm['for_sort']

sns.set(rc = {'figure.figsize':(4.2,2)},font_scale=0.6)
sns.heatmap(dff_toplot_bm.T.loc[cat_order],cmap='PRGn_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/BM_select_pathways.svg')
plt.clf()
plt.close('all')

## Alternate plotting vertically (IGNORE OTHERWISE)

sns.set(rc = {'figure.figsize':(2,4.2)},font_scale=0.6)
sns.heatmap(dff_toplot_bm[cat_order],cmap='PRGn_r',center=0,xticklabels=True)
#plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/BM_select_pathways_alternate.svg')
plt.clf()
plt.close('all')

## for sm vs bm (try and sort for making it more beautiful)

gl_toplot = set(pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/gl_sm_vs_bm_toplot.csv')['GL_pathways'])
sb_toplot = set(pd.read_csv('lymphoid_analysis/csvs/GSEA_shared_pathways/sb_sm_vs_bm_toplot.csv')['SB_pathways'])

comb_toplot = []

for i in list(gl_toplot)+list(sb_toplot):
	if i not in comb_toplot:
		comb_toplot.append(i)
	else:
		pass

## Checking how many of comb-toplot in GL and how many in SB

dff_temp = pd.DataFrame(comb_toplot,columns=['pathway'])
dff_temp['GL_present'] = dff_temp['pathway'].map(lambda x:1 if x in combined_d['gl_sm_vs_gl_bm'].columns.tolist() else 0)
dff_temp['SB_present'] = dff_temp['pathway'].map(lambda x:1 if x in combined_d['sb_sm_vs_sb_bm'].columns.tolist() else 0)

dff_temp=dff_temp.sort_values('GL_present',ascending=False).sort_values('SB_present')

## sort dff_temp based on increased NES within GL, within GL and SB combined and within SB

def nes_extract(gl_grp,sb_grp,pathway):
	if gl_grp==1:
		if sb_grp==0:
			val = combined_d['gl_sm_vs_gl_bm'][pathway].sum().item()
			return val
		else:
			val1 = combined_d['gl_sm_vs_gl_bm'][pathway].sum().item()
			val2 = combined_d['sb_sm_vs_sb_bm'][pathway].sum().item()
			val = (val1 + val2)/2
			return val
	else:
		if sb_grp==1:
			val = combined_d['sb_sm_vs_sb_bm'][pathway].sum().item()
			return val

dff_temp['nes'] = dff_temp.reset_index().apply(lambda x:nes_extract(x['GL_present'],x['SB_present'],x['pathway']),axis=1)

## re-sorting dataframe

dff_temp = pd.concat([dff_temp[dff_temp['GL_present'].isin([1]) & dff_temp['SB_present'].isin([0])].sort_values('nes'),
	dff_temp[dff_temp['GL_present'].isin([1]) & dff_temp['SB_present'].isin([1])].sort_values('nes'),
	dff_temp[dff_temp['GL_present'].isin([0]) & dff_temp['SB_present'].isin([1])].sort_values('nes')],axis=0)

## Make combined dataframe where cell-types have tumor group appended to them
## For pathways not_present, add 0 value for all cells etc. 

gl_notpresent = dff_temp[dff_temp['GL_present'].isin([0])]['pathway'].tolist()
sb_notpresent = dff_temp[dff_temp['SB_present'].isin([0])]['pathway'].tolist()

gl_df1 = combined_d['gl_sm_vs_gl_bm'].copy()
sb_df1 = combined_d['sb_sm_vs_sb_bm'].copy()
gl_df1[gl_notpresent] = pd.DataFrame([[0]*len(gl_notpresent)],index=gl_df1.index)
sb_df1[sb_notpresent] = pd.DataFrame([[0]*len(sb_notpresent)],index=sb_df1.index)

## Now subset dataframe on pathways to be plotted and then plot in correct order

gl_df1['new_index'] = gl_df1.index.map(lambda x:x+'_GL')
sb_df1['new_index'] = sb_df1.index.map(lambda x:x+'_SB')

gl_df1 = gl_df1.set_index('new_index')[dff_temp['pathway']]
sb_df1 = sb_df1.set_index('new_index')[dff_temp['pathway']]

glsb_df = pd.concat([gl_df1,sb_df1],axis=0).copy()

new_order = []
for i in cat_order:
	new_order.append(i+'_SB')
	new_order.append(i+'_GL')

glsb_df.columns = glsb_df.columns.map(lambda x:x.split(' (GO')[0])

sns.set(rc = {'figure.figsize':(7.4,4)},font_scale=0.6)
pol=sns.heatmap(glsb_df.loc[new_order],cmap='PRGn_r',center=0,xticklabels=True)
pol.grid(False)
#plt.tight_layout()
plt.savefig('lymphoid_analysis/gsea_plots/sm_vs_bm/SM_vs_BM_tumors_select_pathways.svg')
plt.clf()
plt.close('all')



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





## Random plotting of Venn Diagrams to show overlap between GL and SB overall across
## cell types

venn_d_sm = {'GL_SM_vs_Sham_SM':set(combined_d['gl_sm_vs_sham_sm'].columns),'SB_SM_vs_Sham_SM':set(combined_d['sb_sm_vs_sham_sm'].columns)}
venn_d_bm = {'GL_BM_vs_Sham_BM':set(combined_d['gl_bm_vs_sham_bm'].columns),'SB_BM_vs_Sham_BM':set(combined_d['sb_bm_vs_sham_bm'].columns)}

from venn import venn

path_to_save = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/gsea_plots/'

fig,ax = plt.subplots(figsize=(4,8),nrows=2)

for i,j,k in zip(['Tumor_SM_vs_Sham_SM','Tumor_BM_vs_Sham_BM'],[venn_d_sm,venn_d_bm],ax.ravel()):
	venn(j, cmap=["deeppink","Green"],ax=k)
	k.title.set_text(i)

plt.tight_layout()
plt.savefig(path_to_save+'combined_overlap_venn_diagram.png',dpi=300)
plt.savefig(path_to_save+'combined_overlap_venn_diagram.svg',dpi=300)
plt.clf()
plt.close('all')

## Now plotting similarly for SM vs BM comparisons

venn_smbm_d = {'Sham_SM_vs_Sham_BM':set(combined_d['sham_sm_vs_sham_bm'].columns),
'GL_SM_vs_GL_BM':set(combined_d['gl_sm_vs_gl_bm'].columns),
'SB_SM_vs_SB_BM':set(combined_d['sb_sm_vs_sb_bm'].columns)}

path_to_save = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/gsea_plots/sm_vs_bm/venn/'

fig,ax = plt.subplots(figsize=(4,5),nrows=1)
venn(venn_smbm_d, cmap=["black","deeppink","Green"],ax=ax)
ax.title.set_text('SM_vs_BM')
plt.tight_layout()
plt.savefig(path_to_save+'combined_overlap_venn_diagram_pathways.png',dpi=300)
plt.savefig(path_to_save+'combined_overlap_venn_diagram_pathways.svg',dpi=300)
plt.clf()
plt.close('all')