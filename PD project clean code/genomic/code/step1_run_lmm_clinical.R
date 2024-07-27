suppressMessages(library('lme4'))
suppressMessages(library('lmerTest'))

df=read.csv('../data/planB_clinical.csv')

traits=c('updrs1','updrs2','updrs3','updrs4','schwab','tremor_scores','pigd_scores','moca','benton',
         'hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall',
         'lns','semantic_fluency','symbol_digit','gds','stai','scopa','rem','ess','quip',
         'gco','alpha_syn','total_tau','abeta_42','p_tau181p')

samples=unique(df$PATNO)

results_intercept <- matrix(nrow = length(traits), ncol = length(samples),dimnames = list(traits, samples))
results_slope <- matrix(nrow = length(traits), ncol = length(samples),dimnames = list(traits, samples))
results_p <- matrix(nrow = length(traits), ncol = 1,dimnames = list(traits))

for (i in 1:length(traits)){
  y=traits[i]
  formular=paste0(y,' ~ age_onset_norm + pd_duration_norm + as.factor(SEX) + as.factor(is_white) + as.factor(is_family_pd) + ledd_norm*time + (1+time|PATNO)')
 if (!all(is.na(df[,paste0(y,'.BL')]))){#check if BL column is all NA (updrs4 is all NA)
    formular=paste0(formular,'+',y,'.BL')
 }
  model=lmer(as.formula(formular),data=df)
  result=ranef(model)$PATNO
  # order results
  ordered_result <- result[match(samples, rownames(result)),]
  # append results
  results_intercept[i, ] <- ordered_result[,1]
  results_slope[i, ] <- ordered_result[,2]
  results_p[i, ] <- summary(model)$coefficients['time','Pr(>|t|)']
  print(paste('finish',y))
}

results_slope_df <- as.data.frame(results_slope)
rownames(results_slope_df) <- traits
colnames(results_slope_df) <- samples

results_intercept_df <- as.data.frame(results_intercept)
rownames(results_intercept_df) <- traits
colnames(results_intercept_df) <- samples

results_p_df <- as.data.frame(results_p)
rownames(results_p_df) <- traits
colnames(results_p_df) <- 'pval'

write.csv(results_slope_df,paste0('../result/planB/trait_random_slope.csv'),row.names = T, quote = F)
write.csv(results_intercept_df,paste0('../result/planB/trait_random_intercept.csv'),row.names = T,quote = F)
write.csv(results_p_df,paste0('../result/planB/trait_pvalue.csv'),row.names = T,quote = F)

