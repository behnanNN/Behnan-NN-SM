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
m2h=pd.read_csv('mouse_to_human_conversion.csv',index_col=0)
m2h_d = dict(zip(m2h['external_gene_name'],m2h['hsapiens_homolog_associated_gene_name']))

## Doing GSEA analysis for all the comparisons
dir_names = ['sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm']

file_list=os.listdir('/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/pydeseq2_analysis/sb_sm_vs_sb_bm/nocutoff')
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

pathways_all.to_csv('csvs/GSEA_analysis/sm_vs_bm/summary_total_differential_pathways_across_groups.csv')
pathways_up.to_csv('csvs/GSEA_analysis/sm_vs_bm/summary_Upregulated_pathways_across_groups.csv')
pathways_dn.to_csv('csvs/GSEA_analysis/sm_vs_bm/summary_Downregulated_pathways_across_groups.csv')

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
	tmpdf.to_csv('csvs/GSEA_analysis/sm_vs_bm/'+'Up_down_pathways_'+i+'.csv')
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

# for i in tqdm(list(adatas_d)):
# 	adata=adatas_d[i]
# 	fig,ax=plt.subplots(figsize=(20,10))
# 	sc.pl.rank_genes_groups_dotplot(adata,show=False,dendrogram=False,ax=ax)
# 	plt.tight_layout()
# 	plt.savefig('pathwys_'+i+'.png',bbox_inches='tight',dpi=300)
# 	plt.clf()
# 	plt.close('all')

## Saving for tomorrow

for i in tqdm(list(adatas_d)):
	adatas_d[i].write('pathways_adata/pathways_'+i+'.h5ad')

## Now plotting venn diagram for each cell type to see what is the overlap of pathways

from venn import venn

cell_names = [x.split('.csv')[0] for x in file_names]

venn_d_total = {}
venn_d_up = {}
venn_d_dn = {}

for i in tqdm(cat_order):
	sham_df = pathways_d['sham_sm_vs_sham_bm'][i].copy()
	gl_df = pathways_d['gl_sm_vs_gl_bm'][i].copy()
	sb_df = pathways_d['sb_sm_vs_sb_bm'][i].copy()
	venn_d_total[i] = {'Sham_Total':set(sham_df['Term']),
	'GL_Total':set(gl_df['Term']),'SB_Total':set(sb_df['Term'])}
	sham_dfup = sham_df[sham_df['NES']>0]
	sham_dfdn = sham_df[sham_df['NES']<0]
	gl_dfup = gl_df[gl_df['NES']>0]
	gl_dfdn = gl_df[gl_df['NES']<0]
	sb_dfup = sb_df[sb_df['NES']>0]
	sb_dfdn = sb_df[sb_df['NES']<0]
	venn_d_up[i] = {'Sham_UP':set(sham_dfup['Term']),
	'GL_UP':set(gl_dfup['Term']),'SB_UP':set(sb_dfup['Term'])}
	venn_d_dn[i] = {'Sham_DN':set(sham_dfdn['Term']),
	'GL_DN':set(gl_dfdn['Term']),'SB_DN':set(sb_dfdn['Term'])}

path_to_save = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/gsea_plots/sm_vs_bm/venn_plots/'

fig,ax = plt.subplots(figsize=(10,18),nrows=5,ncols=3)
fig.subplots_adjust(top=0.8)
for k,axis in zip(list(venn_d_total), ax.ravel()):
	venn(venn_d_total[k], cmap=["Red","Black","Blue"],ax=axis)
	axis.title.set_text(k)

fig.suptitle('Total_Pathways'+'_overlapping_sm_vs_bm', fontsize=16,y=0.99)
plt.tight_layout()
plt.savefig(path_to_save+'Total_overlapping_pathways_venn.png')
plt.savefig(path_to_save+'Total_overlapping_pathways_venn.svg')
plt.clf()
plt.close('all')

fig,ax = plt.subplots(figsize=(10,18),nrows=5,ncols=3)
fig.subplots_adjust(top=0.8)
for k,axis in zip(list(venn_d_up), ax.ravel()):
	venn(venn_d_up[k], cmap=["Red","Black","Blue"],ax=axis)
	axis.title.set_text(k)

fig.suptitle('Upregulated_Pathways'+'_overlapping_sm_vs_bm', fontsize=16,y=0.99)
plt.tight_layout()
plt.savefig(path_to_save+'Upregulated_pathways_overlapping_venn.png')
plt.savefig(path_to_save+'Upregulated_pathways_overlapping_venn.svg')
plt.clf()
plt.close('all')

fig,ax = plt.subplots(figsize=(10,18),nrows=5,ncols=3)
fig.subplots_adjust(top=0.8)
for k,axis in zip(list(venn_d_dn), ax.ravel()):
	venn(venn_d_dn[k], cmap=["Red","Black","Blue"],ax=axis)
	axis.title.set_text(k)

fig.suptitle('Downregulated_Pathways'+'_overlapping_sm_vs_bm', fontsize=16,y=0.99)
plt.tight_layout()
plt.savefig(path_to_save+'Downregulated_overlapping_venn.png')
plt.savefig(path_to_save+'Downregulated_overlapping_venn.svg')
plt.clf()
plt.close('all')



## Starting 12 June 2024 Analysis

## First performing semantic clustering on total pathways on Sham SM vs Sham BM, GL SM vs GL BM,
## and SB SM vs SB BM

dir_names = ['sham_sm_vs_sham_bm','gl_sm_vs_gl_bm','sb_sm_vs_sb_bm']

total_pathways_d = {}

for i in tqdm(dir_names):
	tmpd = pathways_d[i]
	celltypes = list(tmpd)
	tmpls = []
	for j in celltypes:
		df=tmpd[j]
		if df.shape[0]!=0:
			tmpls=tmpls+df['Term'].to_list()
		else:
			pass
	tmpls = list(set(tmpls))
	tmpls = ['GO:' + x.split('(GO:')[1].split(')')[0] for x in tmpls]
	total_pathways_d[i] = pd.DataFrame(tmpls,columns=[i])

total_pathways_d['sham_sm_vs_sham_bm'].to_csv('csvs/For_semantic_clustering/sham_sm_vs_bm_forclustering.csv')
total_pathways_d['gl_sm_vs_gl_bm'].to_csv('csvs/For_semantic_clustering/gl_sm_vs_bm_forclustering.csv')
total_pathways_d['sb_sm_vs_sb_bm'].to_csv('csvs/For_semantic_clustering/sb_sm_vs_bm_forclustering.csv')

## Making a dataframe where all pathways are listed and also listed whether they are common in all or not

all_pathways = [pathways_d[x][y]['Term'].tolist() for x in ['gl_sm_vs_gl_bm', 'sb_sm_vs_sb_bm'] for y in list(pathways_d[x])]
all_pathways = [x for xs in all_pathways for x in xs]
all_pathways = list(set(all_pathways))

new_ls = []

for i in tqdm(all_pathways):
	kol = 0
	tmp_ls = []
	for j in ['gl_sm_vs_gl_bm', 'sb_sm_vs_sb_bm']:
		if i in combined_d[j].columns:
			kol=kol+1
			tmp_ls.append(combined_d[j][i].mean().item())
		else:
			tmp_ls.append(0)
	if kol == 2:
		tmp_ls.append('Shared')
	else:
		tmp_ls.append('Not_Shared')
	new_ls.append(tmp_ls)

all_df = pd.DataFrame(new_ls,columns=['GL_SM_vs_BM','SB_SM_vs_BM','Shared_or_not'],index=all_pathways).sort_values('Shared_or_not',ascending=True)
all_df.to_csv('sm_vs_bm_all_tumor_tolook.csv')


## Looking for common pathways across samples



## Exporting same combined_sm_bm_df with full names rather than GO IDs.

## Importing back semantic clustering results

sem_gl_df = pd.read_csv('GO_semantic_clustering/SM_vs_BM/semantic_cluster_GL_SM_vs_BM.csv')
sem_sb_df = pd.read_csv('GO_semantic_clustering/SM_vs_BM/semantic_cluster_SB_SM_vs_BM.csv')
sem_sham_df = pd.read_csv('GO_semantic_clustering/SM_vs_BM/semantic_cluster_Sham_SM_vs_BM.csv')

sem_sham_d = dict(zip(sem_sham_df['id'],sem_sham_df['cluster']))
sem_gl_d = dict(zip(sem_gl_df['id'],sem_gl_df['cluster']))
sem_sb_d = dict(zip(sem_sb_df['id'],sem_sb_df['cluster']))

new_d = {}

for i in tqdm(dir_names):
	tmpd = pathways_d[i]
	celltypes = list(tmpd)
	tmpls = []
	for j in celltypes:
		df=tmpd[j]
		if df.shape[0]!=0:
			tmpls=tmpls+df['Term'].to_list()
		else:
			pass
	tmpls = list(set(tmpls))
	tmplss = ['GO:' + x.split('(GO:')[1].split(')')[0] for x in tmpls]
	new_d[i] = pd.DataFrame([tmpls,tmplss],index=[i,i+'_Term']).T

new_d['sham_sm_vs_sham_bm']['cluster'] = new_d['sham_sm_vs_sham_bm']['sham_sm_vs_sham_bm_Term'].map(sem_sham_d)
new_d['gl_sm_vs_gl_bm']['cluster'] = new_d['gl_sm_vs_gl_bm']['gl_sm_vs_gl_bm_Term'].map(sem_gl_d)
new_d['sb_sm_vs_sb_bm']['cluster'] = new_d['sb_sm_vs_sb_bm']['sb_sm_vs_sb_bm_Term'].map(sem_sb_d)

new_d['sham_sm_vs_sham_bm'].sort_values('cluster',inplace=True)
new_d['gl_sm_vs_gl_bm'].sort_values('cluster',inplace=True)
new_d['sb_sm_vs_sb_bm'].sort_values('cluster',inplace=True)

new_d['sham_sm_vs_sham_bm'].to_csv('GO_semantic_clustering/SM_vs_BM/sham_processed.csv')
new_d['gl_sm_vs_gl_bm'].to_csv('GO_semantic_clustering/SM_vs_BM/gl_processed.csv')
new_d['sb_sm_vs_sb_bm'].to_csv('GO_semantic_clustering/SM_vs_BM/sb_processed.csv')

## Now plotting selected pathways such that the GL unique are first, 
## then GL and SB shared pathways, and then SB unique pathways. 

gl_toplot = set(pd.read_csv('gl_toplot.csv')['GL_pathways'])
sb_toplot = set(pd.read_csv('sb_toplot.csv')['SB_pathways'])

comb_toplot = []

for i in list(gl_toplot)+list(sb_toplot):
	if i not in comb_toplot:
		comb_toplot.append(i)
	else:
		pass

## Checking how many of comb-toplot in GL and how many in SB

dff_temp = pd.DataFrame(comb_toplot,columns=['pathway'])
dff_temp['GL_present'] = dff_temp['pathway'].map(lambda x:1 if x in new_d['gl_sm_vs_gl_bm']['gl_sm_vs_gl_bm'].tolist() else 0)
dff_temp['SB_present'] = dff_temp['pathway'].map(lambda x:1 if x in new_d['sb_sm_vs_sb_bm']['sb_sm_vs_sb_bm'].tolist() else 0)

dff_temp=dff_temp.sort_values('GL_present',ascending=False).sort_values('SB_present')

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
glsb_df.columns = [x.split(' (GO')[0] for x in glsb_df.columns]

new_order = []
for i in cat_order:
	new_order.append(i+'_SB')
	new_order.append(i+'_GL')


sns.set(rc = {'figure.figsize':(9.4,6)},font_scale=0.6)
pol=sns.heatmap(glsb_df.loc[new_order],cmap='RdBu_r',center=0,xticklabels=True)
pol.grid(False)
#plt.tight_layout()
plt.savefig('gsea_plots/SM_vs_BM/SM_vs_BM_tumors_select_pathways.svg')
plt.clf()
plt.close('all')


## Venn diagram for Sham, GL261, and SB28

from venn import venn

dict_for_venn = {'Sham':set(combined_d['sham_sm_vs_sham_bm'].columns),
'GL261':set(combined_d['gl_sm_vs_gl_bm'].columns),'SB28':set(combined_d['sb_sm_vs_sb_bm'].columns)}

venn(dict_for_venn,cmap=['black','red','green'])
plt.tight_layout()
plt.savefig('gsea_plots/venn_sm_vs_bm.png',dpi=300)
plt.savefig('gsea_plots/venn_sm_vs_bm.svg',dpi=300)
plt.clf()
plt.close('all')