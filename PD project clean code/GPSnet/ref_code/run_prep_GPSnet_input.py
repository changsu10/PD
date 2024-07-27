import os

######### set/change path and parameters

## planB
file_path='../../updrs3_off/result/planB/per_trait_filter100/'
version='matlab'
save_path='../planB/'
VAR_PVAL_NAME='P_Value'
VAR_GENE_NAME='Gene'
VAR_BETA_NAME='Beta'
pval_cutoffs="'[0.1,0.05]'"
adjust_p=False
is_rmOutlier=False
is_norm=True


# ## FMNW_interpolation
# file_path='../result/FM_NW_interpolation_result/'
# save_path='../GPSnet/FMNW_interpolation/'
# VAR_PVAL_NAME='pvalue'
# VAR_GENE_NAME='gene'
# VAR_BETA_NAME='estimate'
# pval_cutoffs="'[0.05,0.1]'"
# adjust_p=False
# is_rmOutlier=False
# is_norm=True
# version='matlab'

cmd='python prep_GPSnet_input.py --file_path=%s --save_path=%s --VAR_PVAL_NAME=%s --VAR_GENE_NAME=%s --VAR_BETA_NAME=%s --pval_cutoffs=%s --version=%s' % \
    (file_path,save_path,VAR_PVAL_NAME,VAR_GENE_NAME,VAR_BETA_NAME,pval_cutoffs,version)
if adjust_p:
    cmd+=' --adjust_p'
if is_rmOutlier:
    cmd += ' --is_rmOutlier'
if is_norm:
    cmd += ' --is_norm'
print(cmd)
os.system(cmd)