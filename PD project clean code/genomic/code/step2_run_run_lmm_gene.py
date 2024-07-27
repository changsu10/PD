import os
import os.path
import sys
import numpy as np
import pandas as pd
import argparse
import subprocess
from multiprocessing import Pool

gene_path='../data/planB_count.csv'
save_path='../result/planB/'

gene_df=pd.read_csv(gene_path,index_col=0)
# pre-filtering as DESeq2
keep = (gene_df >= 10).sum(axis=1) >= 3#keep genes if have at least 3 >10 count
gene_df = gene_df[keep]

num_gene=gene_df.shape[0]
all_sample_names=gene_df.columns.tolist()
all_sample_names=['-'.join(i.split('-')[:2]) for i in all_sample_names]
all_sample_names=list(set(all_sample_names))

####write 100 genes into one file
if not os.path.exists('../result/planB/intermediate_gene/'):
    os.makedirs('../result/planB/intermediate_gene/')

for i in range(num_gene//100+1):
    start=i*100
    end=(i+1)*100
    sub=gene_df.iloc[start:end,:]
    sub.to_csv('../result/planB/intermediate_gene/sub'+str(i)+'.csv')

### 
rcmd='run_lmm_gene.R'
meta_f='../data/planB_meta.csv'
save_f=save_path+'intermediate_ranef/'
if not os.path.exists(save_f):
    os.makedirs(save_f)

n=0
for i in range(num_gene//100+1):
    gene_f='../result/planB/intermediate_gene/sub'+str(i)+'.csv'
    if n%20==0:
        os.system('Rscript %s %s %s %s' % (rcmd,gene_f,meta_f,save_f+'sub'+str(i)))
    else:
        os.system('Rscript %s %s %s %s &' % (rcmd,gene_f,meta_f,save_f+'sub'+str(i)))
    n+=1
    print('finish'+str(i))


slope_df=pd.DataFrame(columns=all_sample_names)
intercept_df=pd.DataFrame(columns=all_sample_names)
p_df=pd.DataFrame(columns=['pval'])
for i in range(num_gene//100+1):
    one_slope=pd.read_csv(save_path+'intermediate_ranef/'+'sub'+str(i)+'_random_slope.csv',index_col=0)
    one_intercept=pd.read_csv(save_path+'intermediate_ranef/'+'sub'+str(i)+'_random_intercept.csv',index_col=0)
    one_p=pd.read_csv(save_path+'intermediate_ranef/'+'sub'+str(i)+'_pvalue.csv',index_col=0)
    # check if missing columns
    if ((one_slope.shape[1]==len(all_sample_names)) and (one_intercept.shape[1]==len(all_sample_names))):
        pass
    else:
        print('unmatched columns in '+str(i))
        slope_missing_cols = set(all_sample_names) - set(one_slope.columns)
        intercept_missing_cols = set(all_sample_names) - set(one_intercept.columns)
        for col in slope_missing_cols:
            one_slope[col] = np.nan
        for col in intercept_missing_cols:
            one_intercept[col] = np.nan
    # reorder
    one_slope = one_slope[all_sample_names]
    one_intercept = one_intercept[all_sample_names]
    # append
    if i == 0:
        slope_df = one_slope
        intercept_df = one_intercept
        p_df = one_p
    else:
        slope_df = pd.concat([slope_df, one_slope], axis=0)
        intercept_df = pd.concat([intercept_df, one_intercept], axis=0)
        p_df = pd.concat([p_df, one_p], axis=0)


slope_df.to_csv(save_path+'gene_random_slope.csv')
intercept_df.to_csv(save_path+'gene_random_intercept.csv')
p_df.to_csv(save_path+'gene_pval.csv')



