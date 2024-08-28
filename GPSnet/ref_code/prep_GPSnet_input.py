##### prepare data used for GPSnet
###### Check if calculate_input() is the way you want
import pandas as pd
import numpy as np
import math
import scipy.io
import copy
from scipy import stats
import os
from statsmodels.stats.multitest import multipletests
from scipy.stats.mstats import winsorize
import argparse
import ast

parser = argparse.ArgumentParser(description='Take inputs for preparing GPSnet input data')
parser.add_argument("--file_path", help="input file path", required=True, type = str)
parser.add_argument("--save_path", help="output file save path", required=True, type = str)

parser.add_argument("--VAR_PVAL_NAME", help="colname in the input file refer to p values", required=True, type = str)
parser.add_argument("--VAR_GENE_NAME", help="colname in the input file refer to gene names", required=True, type = str)
parser.add_argument("--VAR_BETA_NAME", help="colname in the input file refer to beta", required=True, type = str)

parser.add_argument("--pval_cutoffs", help="p value cutoffs", required=True, type = str)

parser.add_argument("--adjust_p", help="whether to adjust p", default=False, action='store_true')
parser.add_argument("--is_rmOutlier", help="whether to remove outliers", default=False, action='store_true')
parser.add_argument("--is_norm", help="whether to normalize input", default=False, action='store_true')

parser.add_argument("--version", help="matlab or python", required=True, type = str)

args = parser.parse_args()

file_path=args.file_path
save_path_parent=args.save_path
VAR_PVAL_NAME=args.VAR_PVAL_NAME
VAR_GENE_NAME=args.VAR_GENE_NAME
VAR_BETA_NAME=args.VAR_BETA_NAME
pval_cutoffs=ast.literal_eval(args.pval_cutoffs)
adjust_p=args.adjust_p
is_rmOutlier=args.is_rmOutlier
is_norm=args.is_norm
version=args.version

################
file_list=os.listdir(file_path)
file_list=[i for i in file_list if i.endswith('.csv')]

suffix=file_path.split('/')[-2]
save_path=save_path_parent+version+'/'+suffix

if adjust_p:
    save_path = save_path + '_padj'
if is_rmOutlier:
    save_path = save_path + '_rmOut'
if is_norm:
    save_path=save_path+'_norm'

save_code_path=save_path+'/code/'
save_path=save_path+'/data/'

if not os.path.exists(save_path):
    os.makedirs(save_path)

print(save_path)

# reference files
gene_id_ref_path='../ref_data/gene_vocab.csv'
gene_id_ref_df=pd.read_csv(gene_id_ref_path)

def calculate_input(one_gene_df):#calculate the input value for each gene
    '''
    input: one_gene_df, sub_df for one gene ncbi id with 3 columns: id, pval, beta
    '''
    if one_gene_df.shape[0]>1:
        print('Multiple matched values for one ncbi id. Return the maximum')
        v_max=0
        for i in range(one_gene_df.shape[0]):
            v=-math.log10(one_gene_df.iloc[i,:]['PVAL']) * abs(one_gene_df.iloc[i,:]['BETA'])
            #v=abs(one_gene_df.iloc[i,:]['BETA'].item())
            if v > v_max:
                v_max=v
    else:
        #v_max = abs(one_gene_df['BETA'].item())  
        v_max=-math.log10(one_gene_df['PVAL']) * abs(one_gene_df['BETA'])   
    if isinstance(v_max, pd.Series):
        v_max=v_max.item()
    return(v_max)

def prep_input(file,p_col,gene_col,beta_col,cal_function,save_path,gene_split=None,is_norm=False,p_cut=0.05,
               adjust_p=False,is_rmOutlier=False,version=version):
    '''
    file: input file path
    p_col: col names for p values (or p adjust)
    gene_col: col names for gene
    gene_split: optional, split character to split genes if multiple ones stored in one line
    beta_col: col names for other input parameters, usually beta or theta
    cal_function: a predefined function to calculate the input (such abs, or -logp)
    is_norm: Ture or False. If true, will use z-score+max-min normalization
    save_path: save path of the output file
    '''
    df=pd.read_csv(file)
    if adjust_p==True:
        p_adj= multipletests(df[p_col], alpha=0.05, method='fdr_bh')[1]
        df[p_col]=p_adj
    ## select significant genes whose p < p_cut and split its genes
    sig_df_list=[]#save gene, p, beta
    for i in range(df.shape[0]):
        one_line=df.iloc[i,]
        pval=one_line[p_col]
        if pval < p_cut:
            gene_str=one_line[gene_col]
            if pd.isnull(gene_str):#no matched gene 
                print('no matched gene for line %i' % i)
            else:
                if gene_split is not None:
                    if gene_split in gene_str:
                        genes=gene_str.split(gene_split)
                    else:#one gene
                        genes=[gene_str]
                else:#one gene
                    genes=[gene_str]
                for g in genes:
                    sig_df_list.append([g,one_line[p_col],one_line[beta_col]])
    sig_df=pd.DataFrame(sig_df_list,columns=['GENE','PVAL','BETA'])
    
    ## map gene name/ensembl ID to NCBI id
    id_df=[]
    not_found_ID=[]
    for i in range(sig_df.shape[0]):
        gene=sig_df.iloc[i,0]
        pval=sig_df.iloc[i,1]
        beta=sig_df.iloc[i,2]
        if gene.startswith('ENSG'):
            if '.' in gene:#if is a transcript
                gene=gene.split('.')[0]
            one_line_ref = gene_id_ref_df.loc[gene_id_ref_df.ensembl_id==gene]
        else:
            one_line_ref=gene_id_ref_df.loc[gene_id_ref_df.symbol==gene]
        if len(one_line_ref)>0:
            ncbi_id=one_line_ref['ncbi_id']
            for j in ncbi_id.values:#if one gene has multiple NCBI ids
                id_df.append([j,pval,beta])
        else:
            not_found_ID.append(gene)

    id_df = pd.DataFrame(id_df, columns=['ncbi_id', 'PVAL', 'BETA'])
    num_gene=id_df.shape[0]
    if num_gene==0:#no sig genes
        return np.nan,0
    
    all_gene_list=list(set(list(id_df.ncbi_id)))
    input_df=[]
    for g in all_gene_list:
        sub_df=id_df.loc[id_df['ncbi_id']==g]
        v = cal_function(sub_df)
        input_df.append([g, v])
    input_df = pd.DataFrame(input_df, columns=['ncbi_id', 'sig'])

    if is_rmOutlier:
        input_df['sig'] = winsorize(input_df['sig'], limits=[0.05, 0.05])

    raw_stats = input_df['sig'].describe()
    # normalize value
    if is_norm:
        norm_value=stats.zscore(input_df['sig'])#z-score
        norm_value=(norm_value-np.min(norm_value))/(np.max(norm_value)-np.min(norm_value))#min-max
        norm_df=copy.deepcopy(input_df)
        norm_df['sig']=norm_value
        save_df=norm_df
    else:
        save_df=input_df

    if version=='matlab':
        scipy.io.savemat(save_path+'.mat', {'Mutation':np.array(save_df)})
        #### prepare gene length mat, set 1 for all genes
        all_gene_id_list = list(set(save_df.ncbi_id))
        gene_length_list = [1] * len(all_gene_id_list)
        gene_length_df = pd.DataFrame({'ncbi_id': all_gene_id_list, 'length': gene_length_list})
        id_mat = np.array(gene_length_df)
        scipy.io.savemat(save_path + '_gene_length_df.mat', {'Gene_Length': id_mat})
        #save_df.to_csv(save_path+'.csv',index=False, header=False)
    else:
        save_df.to_csv(save_path+'.csv',index=False, header=False)

    return raw_stats,num_gene

stats_cols=['trait']
for p_cutoff in pval_cutoffs:
    prefix='pval'+str(p_cutoff)
    stats_cols+=[prefix+'_num',prefix+'_mean',prefix+'_std',prefix+'_min',prefix+'_median',prefix+'_max']

summary_df=pd.DataFrame(columns=stats_cols)

for f in file_list:
    save_name=f.split('.csv')[0].replace(" ", "_").replace("-", "_").replace(".", "_")#replace space and - and . as _
    #has suffix, assume _xx_xx.csv
    new_row_list = [('_').join(save_name.split('_')[:-2])]
    for p_cutoff in pval_cutoffs:
        if not os.path.exists(save_path + 'pval'+str(p_cutoff)+'/'):
            os.makedirs(save_path + 'pval'+str(p_cutoff)+'/')
        raw_stats,num_gene=prep_input(file_path+f,VAR_PVAL_NAME,VAR_GENE_NAME,VAR_BETA_NAME,
                                      calculate_input,save_path+'pval'+str(p_cutoff)+'/'+save_name,
                                      is_norm=is_norm,p_cut=p_cutoff,adjust_p=adjust_p,is_rmOutlier=is_rmOutlier,version=version)
        if num_gene==0:
            new_row_list += [0,np.nan,np.nan,np.nan,np.nan,np.nan]
        else:
            raw_stats = [raw_stats.to_list()[i] for i in [1, 2, 3, 5, 7]]
            new_row_list+=[num_gene]+raw_stats

    new_row_df = pd.DataFrame([new_row_list], columns=stats_cols)
    #summary_df = summary_df.append(new_row_df, ignore_index=True)
    summary_df = pd.concat([summary_df,new_row_df], ignore_index=True)

summary_df.to_csv(save_path+'summary_stats.csv',index=False)

if version=='matlab':
    if not os.path.exists(save_code_path):
        os.makedirs(save_code_path)
    ##### prepare matlab code
    for f in file_list:
        save_name = f.split('.csv')[0].replace(" ", "_").replace("-", "_").replace(".","_")  # replace space and - and . as _
        out = open(save_code_path + 'run_' + save_name + '.m', 'w')
        out.write("trait='%s';\n" % save_name)
        out.write("Raw_Module_Generation(trait,0.5);\n")  # here I use the default hyperparameters
        out.close()
    ##### copy reference code
    copy_cmd=f'cp matlab/* {save_code_path}'
    os.system(copy_cmd)
    ##### copy reference data
    for p_cutoff in pval_cutoffs:
        copy_cmd=f'cp ../ref_data/PPI.mat {save_path}pval{p_cutoff}/'
        os.system(copy_cmd)


