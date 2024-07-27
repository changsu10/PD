## plot score change over time
library(ggplot2)
path='../lmm_result/PD_related/'
df=read.csv('../lmm_result/PD_related/PD_related_combined_data_for_lme.csv')

# m_list=list('rs34637584','rs79987229','rs79987229','rs148534175','rs4266290','rs34778348','rs73038319','rs11949046','rs7134559')
# y_list=list('updrs2','updrs3','pigd_scores','updrs1','p_tau181p','gds','moca','ess','total_tau')

m_list=list('rs2244526','rs8077028','rs17022021','rs8077028','rs115317194','rs149839710','rs186588455')
y_list=list('moca','p_tau181p','updrs2','stai','updrs1','updrs2','schwab')


for (a in 1:length(m_list)){
  #m='rs79987229'
  #y='updrs3'
  m=m_list[a]
  y=y_list[a]
  
  y_index=which(colnames(df)==y)
  m_index=which(colnames(df)==m)
  
  df_sub=df[,c(1:3,y_index,m_index)]
  df_sub$EVENT_ID=as.factor(df_sub$EVENT_ID)
  
  common_events=intersect(as.character(unique(df_sub_w_mutation$EVENT_ID)),as.character(unique(df_sub_wo_mutation$EVENT_ID)))
  common_events=sort(common_events)
  #with mutation
  df_sub_w_mutation=df_sub[df_sub[,5]>0,]
  plot_df=c()
  for (i in common_events){
    v=mean(df_sub_w_mutation[df_sub_w_mutation$EVENT_ID==i,4],na.rm = TRUE)
    plot_df=rbind(plot_df,c(i,v,'w.mutation'))
  }
  
  #without mutation
  df_sub_wo_mutation=df_sub[df_sub[,5]==0,]
  for (i in common_events){
    v=mean(df_sub_wo_mutation[df_sub_wo_mutation$EVENT_ID==i,4],na.rm = TRUE)
    plot_df=rbind(plot_df,c(i,v,'wo.mutation'))
  }
  
  #plot
  plot_df=as.data.frame(plot_df)
  plot_df[plot_df=='NaN']=NA
  plot_df=na.omit(plot_df)
  plot_df[,2]=as.numeric(plot_df[,2])
  colnames(plot_df)=c('EVENT_ID','score','group')
  
  p=ggplot(plot_df, aes(x=EVENT_ID, y=score, colour = group, group = group)) +
    geom_line()+
    ggtitle(paste(y,m,sep='-'))+
    theme_bw()
  
  ggsave(paste0(path,'plot/',y,'-',m,'.png'),p,width=7,height = 5)
  
}
