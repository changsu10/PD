### code to run one model
suppressMessages(library('lme4'))
suppressMessages(library('lmerTest'))
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Need 5 input: 
      trait name, 
      mutation name,
      mutation value,
      trait file,
      save path", call.=FALSE)
}

y=args[1]#trait
y_f=args[2]#trait dataframe
m=args[3]#mutation name
m_f=args[4]#mutation vector
save_path=args[5]#save path

# read data
y_df=read.csv(y_f)
m_df=read.csv(m_f)
df=cbind(y_df,m_df)
df=as.data.frame(df)
colnames(df)=c(colnames(y_df),m)
df$edu_norm=(df$EDUCYRS-mean(df$EDUCYRS))/sd(df$EDUCYRS)

## train model
formular=paste0(y,' ~ edu_norm + age_onset_norm + pd_duration_norm + as.factor(SEX) + as.factor(is_white) + as.factor(is_family_pd) + ',m,'*time + ledd_norm + ledd_norm*time + (1+time|PATNO)')
if (!all(is.na(df[,paste0(y,'.BL')]))){#check if BL column is all NA
  formular=paste0(formular,'+',y,'.BL')
}
model=lmer(as.formula(formular),data=df)


## save m related beta and p
#saveRDS(model, paste0(save_path,'linear_model/',y,'_',m,'.rds'))
pval_m=unname(coef(summary(model))[,'Pr(>|t|)'][m])
pval_mt=unname(coef(summary(model))[,'Pr(>|t|)'][paste0(m,':time')])
beta_m=unname(coef(summary(model))[,'Estimate'][m])
beta_mt=unname(coef(summary(model))[,'Estimate'][paste0(m,':time')])

## 
write.table(pval_m,paste0(save_path,'pval_m/',y,'_',m),row.names = FALSE, col.names = FALSE)
write.table(pval_mt,paste0(save_path,'pval_mt/',y,'_',m),row.names = FALSE, col.names = FALSE)
write.table(beta_m,paste0(save_path,'beta_m/',y,'_',m),row.names = FALSE, col.names = FALSE)
write.table(beta_mt,paste0(save_path,'beta_mt/',y,'_',m),row.names = FALSE, col.names = FALSE)




