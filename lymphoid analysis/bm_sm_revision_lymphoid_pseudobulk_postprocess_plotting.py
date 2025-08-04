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
import seaborn as sns
import os
from tqdm import tqdm
import numpy as np
import random
import pandas as pd
plt.rcParams['svg.fonttype'] = 'none'

## Importing pydeseq2 output files in the form of dictionary

dir_names = os.listdir('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis')
dir_names.remove('summary_plots_overlapping_genes')
file_names = os.listdir('/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/gl_bm_vs_sham_bm/1.2_log2FC_0.05_pval_cutoff')
cell_names = [x.split('.csv')[0] for x in file_names]
common_path = '/gs/gsfs0/users/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/'

master_analyzed_nocutoff = {}
master_analyzed_1_2 = {}

for i in tqdm(dir_names):
	tmp_d = {}
	for j in cell_names:
		tmp_d[j] = pd.read_csv(common_path+i+'/nocutoff/'+j+'.csv',index_col=0)
	master_analyzed_nocutoff[i] = tmp_d	

for i in tqdm(dir_names):
	tmp_d = {}
	for j in cell_names:
		tmp_d[j] = pd.read_csv(common_path+i+'/1.2_log2FC_0.05_pval_cutoff/'+j+'.csv',index_col=0)
	master_analyzed_1_2[i] = tmp_d	

## Let's first plot for each conditions and each cell type, venn diagram depicting upregulated and downregulated genes as well as counts

genes_up_down_d = {}

for i in tqdm(list(master_analyzed_1_2)):
	tmp_d = master_analyzed_1_2[i]
	out_ls = []
	for j in list(tmp_d):
		df=tmp_d[j]
		n_up = df[df['log2FoldChange']>0].shape[0]
		n_dn = df[df['log2FoldChange']<0].shape[0]
		out_ls.append([df.shape[0],n_up,n_dn])
	dff=pd.DataFrame(out_ls,index=list(tmp_d),columns=['Total DEGs 1.2FC','Upregulated','Downregulated'])
	genes_up_down_d[i] = dff

## Plotting this data now

for k in tqdm(list(genes_up_down_d)):
	fig,ax = plt.subplots(figsize=(8,14),nrows=4,ncols=3)
	fig.subplots_adjust(top=0.8)
	for i,j in zip(list(master_analyzed_1_2[k]),ax.ravel()):
		if genes_up_down_d[k].T[i].sum().item() != 0:
			_, _, autotexts=j.pie(genes_up_down_d[k].T[i][['Upregulated','Downregulated']],colors=['Red','Black'],wedgeprops={'alpha':0.7},autopct='%.0f%%')
			j.set_title(i,color='black',fontsize=11,backgroundcolor='sandybrown')
			for autotexts in autotexts:
				autotexts.set_color('white')
		else:
			pass
	fig.legend(['Upregulated Genes','Downregulated Genes'],labelcolor=['red','black'],facecolor='lightgrey',loc='lower center')
	fig.suptitle(k, fontsize=16,y=0.99)
	plt.tight_layout()
	plt.savefig(common_path+'/'+k+'/_summary_up_down_genes.png',dpi=300)
	plt.savefig(common_path+'/'+k+'/_summary_up_down_genes.svg',dpi=300)
	plt.clf()
	plt.close('all')

## Saving Summary data for up and down genes

for i in list(genes_up_down_d):
	genes_up_down_d[i].to_csv('/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/'+i+'/'+i+'_summary_up_down.csv')

## Now looking at commonly UP and commonly Down genes

# Making a venn diagram first

from venn import venn

venn_d = {}

for i,j in tqdm(zip(['gl','sb'],[('gl_bm_vs_sham_bm','gl_sm_vs_sham_sm'),('sb_bm_vs_sham_bm','sb_sm_vs_sham_sm')])):
	bm_d = master_analyzed_1_2[j[0]]
	sm_d = master_analyzed_1_2[j[1]]
	tmp_up_d = {}
	tmp_dn_d = {}
	for k in tqdm(list(bm_d)):		
		bm_df = bm_d[k].copy()
		bm_df.columns = [x+'_bm' for x in bm_df.columns]
		sm_df = sm_d[k].copy()
		sm_df.columns = [x+'_sm' for x in sm_df.columns]
		bm_dfup = bm_df[bm_df['log2FoldChange_bm']>0] 
		bm_dfdn = bm_df[bm_df['log2FoldChange_bm']<0]
		sm_dfup = sm_df[sm_df['log2FoldChange_sm']>0]
		sm_dfdn = sm_df[sm_df['log2FoldChange_sm']<0]
		tmp_up_d[k] = {'bm_up':set(bm_dfup.index),'sm_up':set(sm_dfup.index)}
		tmp_dn_d[k] = {'bm_dn':set(bm_dfdn.index),'sm_dn':set(sm_dfdn.index)}
	venn_d[i] = {'up':tmp_up_d,'dn':tmp_dn_d}

path_to_save = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/summary_plots_overlapping_genes/'

for i in list(venn_d):
	for j in list(venn_d[i]):
		fig,ax = plt.subplots(figsize=(10,18),nrows=4,ncols=3)
		fig.subplots_adjust(top=0.8)
		for k,axis in zip(list(venn_d[i][j]), ax.ravel()):
			venn(venn_d[i][j][k], cmap=["Red","Black"],ax=axis)
			axis.title.set_text(k)
		fig.suptitle(i+'_'+j+'regulated'+'_overlapping_sm_bm', fontsize=16,y=0.99)
		plt.tight_layout()
		plt.savefig(path_to_save+i+'_'+j+'_overlapping_venn.png')
		plt.savefig(path_to_save+i+'_'+j+'_overlapping_venn.svg')
		plt.clf()
		plt.close('all')

venn_d_up = {}
venn_d_dn = {}

for i in tqdm(cell_names):
	sham_df = master_analyzed_1_2['sham_sm_vs_sham_bm'][i].copy()
	sham_df.columns = [x+'_sham' for x in sham_df.columns]
	gl_df = master_analyzed_1_2['gl_sm_vs_gl_bm'][i].copy()
	gl_df.columns = [x+'_gl' for x in gl_df.columns]
	sb_df = master_analyzed_1_2['sb_sm_vs_sb_bm'][i].copy()
	sb_df.columns = [x+'_sb' for x in sb_df.columns]
	sham_dfup = sham_df[sham_df['log2FoldChange_sham']>0] 
	sham_dfdn = sham_df[sham_df['log2FoldChange_sham']<0]
	gl_dfup = gl_df[gl_df['log2FoldChange_gl']>0] 
	gl_dfdn = gl_df[gl_df['log2FoldChange_gl']<0]
	sb_dfup = sb_df[sb_df['log2FoldChange_sb']>0] 
	sb_dfdn = sb_df[sb_df['log2FoldChange_sb']<0]
	venn_d_up[i] = {'Sham_UP':set(sham_dfup.index),'GL_UP':set(gl_dfup.index),'SB_UP':set(sb_dfup.index)}
	venn_d_dn[i] = {'Sham_DN':set(sham_dfdn.index),'GL_DN':set(gl_dfdn.index),'SB_DN':set(sb_dfdn.index)}

path_to_save = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_analysis_new/lymphoid_analysis/pydeseq2_analysis/summary_plots_overlapping_genes/sm_vs_bm/'

fig,ax = plt.subplots(figsize=(10,18),nrows=4,ncols=3)
fig.subplots_adjust(top=0.8)
for k,axis in zip(list(venn_d_up), ax.ravel()):
	venn(venn_d_up[k], cmap=["Red","Black","Blue"],ax=axis)
	axis.title.set_text(k)

fig.suptitle('upregulated'+'_overlapping_sm_vs_bm', fontsize=16,y=0.99)
plt.tight_layout()
plt.savefig(path_to_save+'Upregulated_overlapping_venn.png')
plt.savefig(path_to_save+'Upregulated_overlapping_venn.svg')
plt.clf()
plt.close('all')

fig,ax = plt.subplots(figsize=(10,18),nrows=4,ncols=3)
fig.subplots_adjust(top=0.8)
for k,axis in zip(list(venn_d_dn), ax.ravel()):
	venn(venn_d_dn[k], cmap=["Red","Black","Blue"],ax=axis)
	axis.title.set_text(k)

fig.suptitle('upregulated'+'_overlapping_sm_vs_bm', fontsize=16,y=0.99)
plt.tight_layout()
plt.savefig(path_to_save+'Downregulated_overlapping_venn.png')
plt.savefig(path_to_save+'Downregulated_overlapping_venn.svg')
plt.clf()
plt.close('all')


master_overlap = {}
summary_overlap_ = {}

for i,j in zip(['gl','sb'],[('gl_bm_vs_sham_bm','gl_sm_vs_sham_sm'),('sb_bm_vs_sham_bm','sb_sm_vs_sham_sm')]):
	bm_d = master_analyzed_1_2[j[0]]
	sm_d = master_analyzed_1_2[j[1]]
	tmp_d = {}
	tmp_ls = []
	for k in tqdm(list(bm_d)):
		bm_df = bm_d[k]
		sm_df = sm_d[k]
		bm_dfup = bm_df[bm_df['log2FoldChange_bm']>0] 
		bm_dfdn = bm_df[bm_df['log2FoldChange_bm']<0]
		sm_dfup = sm_df[sm_df['log2FoldChange_sm']>0]
		sm_dfdn = sm_df[sm_df['log2FoldChange_sm']<0]
		common_up_ls = list(set(bm_dfup.index).intersection(set(sm_dfup.index)))
		common_dn_ls = list(set(bm_dfdn.index).intersection(set(sm_dfdn.index)))
		common_df_up = pd.concat([bm_dfup.loc[common_up_ls],sm_dfup.loc[common_up_ls]],axis=0)
		common_df_dn = pd.concat([bm_dfdn.loc[common_up_dn],sm_dfdn.loc[common_dn_ls]],axis=0)
		disconcordant_ls = list(set(bm_df.index).union(set(sm_df.index)))-(set(common_up_ls).union(set(common_dn_ls)))
		disconcordant_df = pd.concat([bm_df.loc[disconcordant_ls],sm_df.loc[disconcordant_ls]],axis=0)
		tmp_d['common_up'] = common_df_up
		tmp_d['common_dn'] = common_df_dn
		tmp_d['disconcordant'] = disconcordant_df
		tmp_ls.append([bm_df.shape[0]])
	master_overlap[i] = tmp_d

## Random plotting of Venn Diagrams to show overlap between GL and SB overall across
## cell types

combined_genes_d = {}

for i in list(master_analyzed_1_2):
	tmp_d = master_analyzed_1_2[i]
	all_genes = [x for k in list(tmp_d) for x in tmp_d[k].index.unique().tolist()]
	all_genes = list(set(all_genes))
	tmpls = []
	for cell_types in list(tmp_d):
		tmpdf = tmp_d[cell_types]
		log2fc_d = dict(zip(tmp_d[cell_types].index,tmp_d[cell_types]['log2FoldChange']))
		log2fc_d.update(dict.fromkeys(list(set(all_genes)-set(list(log2fc_d))),0))
		tmpls.append([log2fc_d[x] for x in all_genes])
	outdf=pd.DataFrame(tmpls,index=list(tmp_d),columns=all_genes)
	combined_genes_d[i] = outdf


## Venn diagram for Tumor SM vs Sham SM

from venn import venn

sm_gl_vs_sham = [master_analyzed_1_2['gl_sm_vs_sham_sm'][x].index.tolist() for x in list(master_analyzed_1_2['gl_sm_vs_sham_sm'])]
sm_gl_vs_sham = set([x for xs in sm_gl_vs_sham for x in xs])
sm_sb_vs_sham = [master_analyzed_1_2['sb_sm_vs_sham_sm'][x].index.tolist() for x in list(master_analyzed_1_2['sb_sm_vs_sham_sm'])]
sm_sb_vs_sham = set([x for xs in sm_sb_vs_sham for x in xs])

dict_for_venn = {'GL261':sm_gl_vs_sham,
'SB28':sm_sb_vs_sham}

venn(dict_for_venn,cmap=['red','green'])
plt.tight_layout()
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_tumor_sm_vs_sham_sm.png',dpi=300)
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_tumor_sm_vs_sham_sm.svg',dpi=300)
plt.clf()
plt.close('all')

## Venn diagram for Tumor BM vs Sham BM

from venn import venn

bm_gl_vs_sham = [master_analyzed_1_2['gl_bm_vs_sham_bm'][x].index.tolist() for x in list(master_analyzed_1_2['gl_bm_vs_sham_bm'])]
bm_gl_vs_sham = set([x for xs in bm_gl_vs_sham for x in xs])
bm_sb_vs_sham = [master_analyzed_1_2['sb_bm_vs_sham_bm'][x].index.tolist() for x in list(master_analyzed_1_2['sb_bm_vs_sham_bm'])]
bm_sb_vs_sham = set([x for xs in bm_sb_vs_sham for x in xs])

dict_for_venn = {'GL261':bm_gl_vs_sham,
'SB28':bm_sb_vs_sham}

venn(dict_for_venn,cmap=['red','green'])
plt.tight_layout()
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_tumor_bm_vs_sham_bm.png',dpi=300)
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_tumor_bm_vs_sham_bm.svg',dpi=300)
plt.clf()
plt.close('all')

## Venn diagram for Sham, GL261 and SB28

sham_sm_vs_bm = [master_analyzed_1_2['sham_sm_vs_sham_bm'][x].index.tolist() for x in list(master_analyzed_1_2['sham_sm_vs_sham_bm'])]
sham_sm_vs_bm = set([x for xs in sham_sm_vs_bm for x in xs])

gl_sm_vs_bm = [master_analyzed_1_2['gl_sm_vs_gl_bm'][x].index.tolist() for x in list(master_analyzed_1_2['gl_sm_vs_gl_bm'])]
gl_sm_vs_bm = set([x for xs in gl_sm_vs_bm for x in xs])

sb_sm_vs_bm = [master_analyzed_1_2['sb_sm_vs_sb_bm'][x].index.tolist() for x in list(master_analyzed_1_2['sb_sm_vs_sb_bm'])]
sb_sm_vs_bm = set([x for xs in sb_sm_vs_bm for x in xs])

dict_for_venn = {'Sham':sham_sm_vs_bm,
'GL261':gl_sm_vs_bm,'SB28':sb_sm_vs_bm}

venn(dict_for_venn,cmap=['black','red','green'])
plt.tight_layout()
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_sm_vs_bm.png',dpi=300)
plt.savefig('lymphoid_analysis/pydeseq2_analysis/venn_sm_vs_bm.svg',dpi=300)
plt.clf()
plt.close('all')