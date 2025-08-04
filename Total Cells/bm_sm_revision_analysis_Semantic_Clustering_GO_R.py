## Performing Semantic clustering using SimplifyEnrichment Package

## Moving to R now

setwd("Z:/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/csvs/For_semantic_clustering")

library(simplifyEnrichment)

goids_df = read.csv("Shared_pathways_for_clustering.csv")

smup_df = read.csv('sm_up_forclustering.csv')
smdn_df = read.csv('sm_dn_forclustering.csv')
bmup_df = read.csv('bm_up_forclustering.csv')
bmdn_df = read.csv('bm_dn_forclustering.csv')

smup_mat = GO_similarity(smup_df$sm_up)
smdn_mat = GO_similarity(smdn_df$sm_dn)
bmup_mat = GO_similarity(bmup_df$bm_up)
bmdn_mat = GO_similarity(bmdn_df$bm_dn)

smup_df = simplifyGO(smup_mat)
write.csv(smup_df,'semantic_cluster_SM_UP.csv')

> table(smup_df$cluster)

  1   2   3   4   5   6   7   8   9  10  11 
 72 109  25   4   1  19   6  28   4   1   1

smdn_df = simplifyGO(smdn_mat)
write.csv(smdn_df,'semantic_cluster_SM_DN.csv')

> table(smdn_df$cluster)

 1  2  3  4  5  6  7  8  9 10 11 12 
19  7 23  5 34  3  2  7  1  5  1  1

bmup_df = simplifyGO(bmup_mat)
write.csv(bmup_df,'semantic_cluster_BM_UP.csv')

> table(bmup_df$cluster)

 1  2  3  4  5  6  7  8  9 10 11 
11  1  6  7  6  1  1  1  1  2  1

bmdn_df = simplifyGO(bmdn_mat)
write.csv(bmdn_df,'semantic_cluster_BM_DN.csv')

> table(bmdn_df$cluster)

  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
149  84  10  74  16   3   1   1   8  10   7   1  15   1   1   1 

## Going back to python and putting Pathway Name in front of the GO_IDs.

import gseapy as gp

lib = gp.get_library('GO_Biological_Process_2023')
lib_ls = list(lib)

go_d = {}

for i in lib_ls:
	go_d[i.split('(')[1].split(')')[0]] = i

## Now importing the files, updating the GO-Term and saving them back

path_to_open = '/gs/gsfs0/home/abdubey/pw/storage/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/csvs/For_semantic_clustering/' 

for i in ['semantic_cluster_BM_DN.csv', 'semantic_cluster_BM_UP.csv', 'semantic_cluster_SM_DN.csv', 'semantic_cluster_SM_UP.csv']:
	dff = pd.read_csv(path_to_open+i,index_col=0)
	dff['Term'] = dff['id'].map(go_d)
	dff[['id','Term','cluster']].to_csv(path_to_open+i)


## Semantic Clustering on SM vs BM comparisons

setwd("Z:/scRNAseq/2024_Analysis/scanpy_analysis/bm_sm_revision_analysis/csvs/For_semantic_clustering")

library(simplifyEnrichment)

sham_smbm = read.csv("sham_sm_vs_bm_forclustering.csv")
gl_smbm = read.csv("gl_sm_vs_bm_forclustering.csv")
sb_smbm = read.csv("sb_sm_vs_bm_forclustering.csv")

sham_mat = GO_similarity(sham_smbm$sham_sm_vs_sham_bm)
gl_mat = GO_similarity(gl_smbm$gl_sm_vs_gl_bm)
sb_mat = GO_similarity(sb_smbm$sb_sm_vs_sb_bm)

sham_df = simplifyGO(sham_mat)
write.csv(sham_df,'semantic_cluster_Sham_SM_vs_BM.csv')

> table(sham_df$cluster)

  1   2   3   4   5   6   7   8   9  10 
288 204  78 166  81   8   8  37   5  21 

gl_df = simplifyGO(gl_mat)
write.csv(gl_df,'semantic_cluster_GL_SM_vs_BM.csv')

> table(gl_df$cluster)

  1   2   3   4   5   6   7 
212 178  90  83  39   5  28

sb_df = simplifyGO(sb_mat)
write.csv(sb_df,'semantic_cluster_SB_SM_vs_BM.csv')

> table(sb_df$cluster)

  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
138 118  68 153  71  55   6  17   7   3   5   2   7   3   6   2   1   1   2 

## Lymphoid new one

gl_mat = GO_similarity(gldf$go_ids)
sb_mat = GO_similarity(sbdf$go_ids)

gl_df = simplifyGO(gl_mat)
sb_df = simplifyGO(sb_mat)

write.csv(gl_df,'semantic_cluster_GL_SM_vs_BM.csv')

> table(sb_df$cluster)
1  2  3  4  5  6  7  8  9 10 11 12 
15 26 37  9  7 15  6  1  2  2  1  1
