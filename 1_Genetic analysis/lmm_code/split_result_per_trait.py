import pandas as pd
import statsmodels.stats.multitest as smm
import numpy as np
from pandas.util.testing import assert_frame_equal
import os
######### Load data
df_pd=pd.read_csv('../lmm_results/PD_combined/linear_pval_mt_combined_all.csv',index_col=0)
beta_pd=pd.read_csv('../lmm_results/PD_combined/linear_beta_mt_combined_all.csv',index_col=0)
map_ref=pd.read_csv('../../../ref/combined_map_gene2snp.csv')

save_path='../lmm_results/PD_combined/per_trait/'
if not os.path.exists(save_path):
        os.makedirs(save_path)
######### user provided trait list
# traits_nonCognition_list=['updrs2','updrs3','schwab',#motor
#                           'updrs1','updrs4','tremor_scores','pigd_scores',# general
#                           'gds','stai',#mood
#                           'scopa',#automatic
#                           'rem','ess',#sleep
#                           'quip','gco',
#                           'alpha_syn','total_tau','abeta_42','p_tau181p'#biomarker
#                          ]
# traits_cognition_list=['moca','benton','hvlt_delayed_recall','hvlt_recog_disc_index',
#                        'hvlt_retention','hvlt_total_recall','lns','semantic_fluency','symbol_digit']
# traits_list=traits_nonCognition_list+traits_cognition_list
traits_list=['updrs3']

#####
# take_nonNA = lambda s1, s2: s2 if s1 is np.nan else s1

###### Save results per trait
for trait in traits_list:
    p_mt=df_pd[trait].to_frame()
    beta_mt=beta_pd[trait].to_frame()
    
    # select the non-NA 
    p_mt=p_mt[p_mt[trait].notna()]#remove NA
    p_mt.rename(columns={trait: 'p_mt'}, inplace=True)
    
    beta_mt=beta_mt[beta_mt[trait].notna()]#remove NA
    beta_mt.rename(columns={trait: 'beta_mt'}, inplace=True)
    
    df_mt = pd.concat([beta_mt, p_mt], axis=1)
    
    #correct pvalue
    p_adj=smm.multipletests(list(df_mt['p_mt']), method='fdr_bh')[1]
    df_mt['p.adj_mt']=p_adj
    df_mt=df_mt.sort_values(['p_mt'])#sort by p value

    l=list(df_mt.index)#SNP list
    #check if p_mt has dulicated SNPs
    if len(l)!=len(set(l)):
        print('!!ERROR!!%s has duplicated SNPs' % trait)

    #map SNP to genes
    gene_sub=[]
    for s in l:
        snp=map_ref.loc[map_ref.SNP==s]
        genes=list(snp.GENES)
        if len(genes)!=len(set(genes)):
            print('trait %s has duplicated genes' % trait)#no duplicated genes
        gene_sub.append([s,'|'.join(genes)])
    gene_sub=pd.DataFrame(gene_sub,columns=['SNP','Genes'])
    gene_sub.set_index('SNP', inplace=True)
    
    #check gene_sub and previous p_mt dimension (if some SNPs don't have genes)
    if gene_sub.shape[0]!=p_mt.shape[0]:
        print('In %s, p_mt and gene_sub dim different' % trait)
    
    # save SNP, beta, p, p.adj, gene
    df_save = pd.concat([df_mt, gene_sub], axis=1)#concat by index
    df_save=df_save.rename_axis('SNP').reset_index()
    
    df_save.to_csv(save_path+trait+'_linear_mt.csv',index=False)
