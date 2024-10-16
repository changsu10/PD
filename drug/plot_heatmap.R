library(data.table)
library(ggplot2)
library(dplyr)
setwd('~/Desktop/amp_pd/drug_repurposing')

score_cutoff=-2
p_cutoff=0.05
save_plot_path='./plot/v2_-2_0.05/'

all_files=list.files('./data/')
combined_df=c()
combined_df_all=c()
for (file in all_files){
  df=fread(paste0('./data/',file))
  df=na.omit(df)
  if (dim(df)[1]>0){
    colnames(df)=c('drug','score','pval')
    sig_df=df[(df$score< score_cutoff) & (df$pval< p_cutoff),]
    trait_name=sub("240721_(.*)_network_proximity.txt", "\\1", file)
    if (trait_name=='motor_and_cognition1_overlapped'){
      trait_name='motor_cog1_common'
    }
    if (trait_name=='motor_and_cognition1'){
      trait_name='motor_cog1_union'
    }
    sig_df$trait=rep(trait_name,dim(sig_df)[1])
    combined_df=rbind(combined_df,sig_df)
    
    df$trait=rep(trait_name,dim(df)[1])
    combined_df_all=rbind(combined_df_all,df)
  }
}

#---- Heatmap -- row:drug, col: traits, color:score, size:p ----
combined_df <- combined_df %>%
  arrange(score, pval)
combined_df$trait <- factor(combined_df$trait, 
                            levels = c('motor','cognition1','cognition2',"motor_cog1_union","motor_cog1_common",
                                       'mood','sleep','autonomic',
                                       'abeta_42','alpha_syn','p_tau181p','total_tau','all_combined'))

# Plot 1
p=ggplot(combined_df, aes(trait, drug, color = score, size = -log10(pval))) +
  geom_point() +
  scale_size_continuous(range = c(2,4))+
  theme_bw() +
  scale_colour_gradient(low = "#132B43",high = "#56B1F7",)+
  theme(legend.position = "right", 
        #panel.grid.major = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(paste0(save_plot_path,'heatmap_all.pdf'),p,width=10,heigh=40,limitsize = FALSE)

# Plot 2
sub_df=combined_df[(combined_df$trait %in% c('motor','cognition1','cognition2',"motor_cog1_union","motor_cog1_common")),]
p=ggplot(sub_df, aes(trait, drug, color = score, size = -log10(pval))) +
  geom_point() +
  scale_size_continuous(range = c(2,4))+
  theme_bw() +
  scale_colour_gradient(low = "#132B43",high = "#56B1F7",)+
  theme(legend.position = "right", 
        #panel.grid.major = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(paste0(save_plot_path,'heatmap_motor_cog.pdf'),p,width=10,heigh=40,limitsize = FALSE)

# Plot 3
sub_df=combined_df[(combined_df$trait %in% c('motor','cognition1','cognition2',"mood","sleep",'autonomic','abeta_42','alpha_syn','total_tau','p_tau181p')),]
p=ggplot(sub_df, aes(trait, drug, color = score, size = -log10(pval))) +
  geom_point() +
  theme_bw() +
  scale_colour_gradient(low = "#132B43",high = "#56B1F7",)+
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(paste0(save_plot_path,'heatmap_sub.pdf'),p,width=10,heigh=50,limitsize = FALSE)

#---- Extract unique (in 1 trait) and common (in >5 traits) drugs ----
combined_df_sub_traits = combined_df[(combined_df$trait %in% c('motor','cognition1','cognition2','mood','sleep','autonomic')),]
count_df <- combined_df_sub_traits %>%
  group_by(drug) %>%
  summarize(trait_count = n_distinct(trait), 
            traits = paste(unique(trait), collapse = "|"))
write.table(count_df,paste0(save_plot_path,'drug_in_noBioMtrait_count.tsv'),sep='\t',row.names = F,quote = F)


#---- Plot selected drug ----
sub_df=combined_df_all[(combined_df_all$trait %in% c('motor','cognition1','cognition2','mood','sleep','autonomic')),]
selected_drug=read.csv('selected_drug.csv',header=F)
sub_df=sub_df[sub_df$drug %in% selected_drug$V1,]

sub_df$trait[sub_df$trait=='cognition2']='moca'
sub_df$trait[sub_df$trait=='cognition1']='cognition'
sub_df$trait <- factor(sub_df$trait, levels = c('motor','cognition','moca','mood','sleep','autonomic'))
sub_df$drug = factor(sub_df$drug,levels=rev(selected_drug$V1))
sub_df$z_score=sub_df$score
p=ggplot(sub_df, aes(trait, drug, color = z_score, size = -log10(pval))) +
  geom_point() +
  scale_size_continuous(range = c(1,7))+
  theme_bw() +
  #scale_color_gradientn(breaks=c(-Inf, -3,0,3, Inf), colors=c("blue","white",'red',"red","red"))+
  scale_color_gradient2(low='#2983b1',mid = "white",high = "#D53389",midpoint = 0,limits=c(-3,3),na.value = "#2e59a7")+#lower the better
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave('heatmap_sub_selected_drug.png',p,width=7,heigh=13.5)



