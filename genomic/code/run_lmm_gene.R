### code to run one model
# input: 100 genes data; meta data

suppressMessages(library('lme4'))
suppressMessages(library('lmerTest'))
suppressMessages(library('dplyr'))
suppressMessages(library('tidyr'))

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Need 3 input: 
      gene file path, 
      meta file path,
      save path", call.=FALSE)
}

g_f=args[1]#
meta_f=args[2]#meta file
save_path=args[3]#save path

# read data
y_df=read.csv(g_f, check.names = FALSE, row.names = 1)
genes=rownames(y_df)

x_df=read.csv(meta_f, row.names = 1)
x_df$Sex=factor(x_df$Sex)
## change race to binary: white vs non-white
x_df$Race <- ifelse(x_df$Race == "White", "white", "non-white")
x_df$Race=factor(x_df$Race)
x_df$time=x_df$Timepoint/12

## make meta the same order as count
x_df <- x_df[match(colnames(y_df), rownames(x_df)), ]
df=cbind(x_df,t(y_df))
df=as.data.frame(df)
#df=df[df$Group=='PD',]
samples=unique(df$PID)

results_intercept <- matrix(nrow = length(genes), ncol = length(samples),dimnames = list(genes, samples))
results_slope <- matrix(nrow = length(genes), ncol = length(samples),dimnames = list(genes, samples))
results_p <- matrix(nrow = length(genes), ncol = 1,dimnames = list(genes))
for (i in 1:length(genes)){
  g=genes[i]
  formular=paste0(g,' ~ age_onset_norm + pd_duration_norm+ Sex + Race  + is_family_pd + ledd_norm*time + (1+time|PID)')
  model=lmer(as.formula(formular),data=df)
  # get random effect
  result=ranef(model)$PID
  # order results
  ordered_result <- result[match(samples, rownames(result)),]
  # append results
  results_intercept[i, ] <- ordered_result[,1]
  results_slope[i, ] <- ordered_result[,2]
  results_p[i, ] <- summary(model)$coefficients['time','Pr(>|t|)']
}

results_slope_df <- as.data.frame(results_slope)
rownames(results_slope_df) <- genes
colnames(results_slope_df) <- samples

results_intercept_df <- as.data.frame(results_intercept)
rownames(results_intercept_df) <- genes
colnames(results_intercept_df) <- samples

results_p_df <- as.data.frame(results_p)
rownames(results_p_df) <- genes
colnames(results_p_df) <- 'pval'

write.csv(results_slope_df,paste0(save_path,'_random_slope.csv'),row.names = T, quote = F)
write.csv(results_intercept_df,paste0(save_path,'_random_intercept.csv'),row.names = T,quote = F)
write.csv(results_p_df,paste0(save_path,'_pvalue.csv'),row.names = T,quote = F)


