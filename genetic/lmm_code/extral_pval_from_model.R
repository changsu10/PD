
f='combined_data_for_lme.csv'
df=read.csv(f)
col_names=colnames(df)

traits=c('updrs2','updrs3','updrs4','schwab','tremor_scores','pigd_scores','moca','benton',
        'hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall',
        'lns','semantic_fluency','symbol_digit','gds','stai','scopa','rem','ess','quip',
        'updrs1','gco','alpha_syn','total_tau','abeta_42','p_tau181p')

mutations=c()
for (i in col_names){
    if (startsWith(i,'rs') | startsWith(i,'chr')){
        mutations=c(mutations,i)
    }
}

pval_all=c()
for (y in traits){
    y_pval=c()
    for (m in mutations){
        mode_name=paste0('./model/',y,'_',m,'.rds')
        if (file.exists(mode_name)){
             model=readRDS(mode_name)
            pval=unname(coef(summary(model))[,'Pr(>|t|)'][m])
            y_pval=c(y_pval,pval)
        } else{
            print(paste('no model for',y,m))
            y_pval=c(y_pval,NA)
        }  
    }
    pval_all=cbind(pval_all,y_pval)
}

pval_all=as.data.frame(pval_all)
colnames(pval_all)=traits
rownames(pval_all)=mutations

write.csv(pval_all,'./pval/all_pval_linear_m.csv',quote=FALSE)



