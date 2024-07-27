import os
import os.path
import sys
import numpy as np
import pandas as pd
import argparse
import subprocess
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument("--snp_path", type=str, default='./PD_related/PD_related_combined_GWAS_for_lme.csv', help='Path of SNP data.')
parser.add_argument("--model_type", type=str,default='linear',help='linear')
parser.add_argument("--snp_type", type=str,default='PD_related',help='PD or PD_related or PD_combined')

args = parser.parse_args()

model_type=args.model_type
snp_path=args.snp_path
snp_type=args.snp_type

traits_cognition_list=['moca','benton','hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall','lns','semantic_fluency','symbol_digit']
###### user provided input traits list
# traits=['updrs2','updrs3','updrs4','schwab','tremor_scores','pigd_scores','moca','benton',
#          'hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall',
#          'lns','semantic_fluency','symbol_digit','gds','stai','scopa','rem','ess','quip',
#          'updrs1','gco','alpha_syn','total_tau','abeta_42','p_tau181p']
traits=['updrs3']#updrs3 only

mutation_df=pd.read_csv(snp_path)
mutations=list(mutation_df.columns)

####write one snp into one file
if not os.path.exists('./'+snp_type+'/GWAS/'):
    os.makedirs('./'+snp_type+'/GWAS/')
    for m in mutations:
        one_m=mutation_df[m]
        one_m.to_csv('./'+snp_type+'/GWAS/'+m+'.csv',index=False)

### Rscript name, save path name
if model_type=='linear':
    save_path='./'+snp_type+'/'+model_type+'_'
else:
    print('!!ERROR!!WRONG model_type')

### create save path sub folders
save_path_types=[save_path+'pval_m',save_path+'pval_mt',save_path+'beta_m',save_path+'beta_mt']
if model_type=='nonlinear':
    save_path_types+=[save_path+'pval_mt2',save_path+'beta_mt2']

for s in save_path_types:
    if not os.path.exists(s):
        os.makedirs(s)

### clinical file name
trait_f='../preprocessed_data/'+snp_type+'_combined_nonGWAS_for_lme.csv'

###
n=1
for t in traits:
    # if trait is congition, add EDU variable
    if t in traits_cognition_list:
        rcmd='one_lmm_cognition_add_edu.R'
    else:
        rcmd='one_lmm.R'
    # loop for each mutations
    for m in mutations:
        m_f='./'+snp_type+'/GWAS/'+m+'.csv'

        if n%50==0:
            os.system('Rscript %s %s %s %s %s %s' % (rcmd,t,trait_f,m,m_f,save_path))
        else:
            os.system('Rscript %s %s %s %s %s %s &' % (rcmd,t,trait_f,m,m_f,save_path))
        n+=1
    
    print('finish'+t)


#### combine estimates, p values
if not os.path.exists('../lmm_result/'+snp_type+'/'):
    os.makedirs('../lmm_result/'+snp_type+'/')

for s in save_path_types:
    df=pd.DataFrame(index=mutations,columns=traits)

    for t in traits:
        for m in mutations:
            v_f_name=s+'/'+t+'_'+m
            if os.path.exists(v_f_name):
                v_f=open(v_f_name,'r')
                v=v_f.readlines()[0].rstrip()
                if v=='NA':
                    df.at[m, t] = np.nan
                else:
                    v=float(v)
                    df.at[m, t] = v
            else:
                df.at[m, t] = np.nan
    name=s.split('/')[-1]
    save_name='../lmm_result/'+snp_type+'/'+name+'_combined_all.csv'
    df.to_csv(save_name)
    #os.system('rm -r %s' % s)

#os.system('rm -r %s' % ('./'+snp_type+'/GWAS/'))













