### save results per trait
#setwd("/Users/manage/Desktop/amp_pd/rnaseq/code")

beta_f='../result/planB/beta_result.csv'
pval_f='../result/planB/pvalue_result.csv'
beta_df=read.csv(beta_f, check.names = FALSE, row.names = 1)
pval_df=read.csv(pval_f, check.names = FALSE, row.names = 1)
save_path='../result/planB/per_trait_filter100/'

if (!dir.exists(save_path)) {
  dir.create(save_path)
}

save_per_trait=function(sig_gene,beta_df,pval_df,save_path,suffix='_beta_pval.csv'){
  traits=colnames(beta_df)
  beta_df=beta_df[sig_gene,, drop = FALSE]
  pval_df=pval_df[sig_gene,, drop = FALSE]
  summary_df=data.frame(row.names=traits)
  sig_num=c()
  for (t in traits){
    beta=beta_df[,t]
    pval=pval_df[,t]
    df=data.frame('Gene'=sig_gene,'Beta'=beta,'P_Value'=pval)
    write.csv(df,paste0(save_path,t,suffix),
              quote=FALSE,row.names =FALSE)
    sig_num=c(sig_num,dim(df[df$P_Value<0.05,])[1])
  }
  summary_df=data.frame('Trait'=traits,'sigNum'=sig_num)
  return(summary_df)
}

################### filter based on raw count data
count_f='../data/planB_count.csv'
count_df=read.csv(count_f, check.names = FALSE, row.names = 1)
sig_gene1=rownames(count_df)[rowSums(count_df>=100)>=500]

sig_gene=sig_gene1
df1=save_per_trait(sig_gene,beta_df,pval_df,
               save_path=save_path,suffix='_beta_pval.csv')#count filter
write.csv(df1,'../result/planB/per_trait_filter100_summary_stats.csv',quote=FALSE,row.names =FALSE)

# ################### filter based on gene associated pvalues
# gene_pval_f='../result/planB/gene_pval.csv'#use for filtering
# deseq2_f="../result/deseq2_sig_genes.csv"

# gene_pval_df=read.csv(gene_pval_f, check.names = FALSE, row.names = 1)
# deseq2_df=read.csv(deseq2_f, check.names = FALSE, row.names = 1)

# # keep genes with p<0.1
# sig_gene=rownames(gene_pval_df)[gene_pval_df<0.1]
# df1=save_per_trait(sig_gene,beta_df,pval_df,
#                    save_path='../result/planB/per_trait_p/',suffix='_beta_pval.csv')

