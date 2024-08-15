library(circlize)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Need 3 input: 
      read path, 
      save path,
      traits", call.=FALSE)
}

file_path=args[1]
save_path=args[2]
traits <- as.vector(strsplit(args[3], ",")[[1]])

top_num=30
for (trait in traits){
  df=read.table(paste0(file_path,trait,'_enrichedPathway_filtered.tsv'),sep='\t',header=T)
  df=df[,c('term_name','p_value','source','term_size')]
  df <- df %>%
    mutate(term_name = ifelse(nchar(term_name) > 80, substr(term_name, 1, 80), term_name))
  ####### bubble plot
  plot_df=df
  plot_df=plot_df[order(plot_df$p_value),][1:top_num,]
  plot_df$logp=-log10(plot_df$p_value)
  
  plot_df=plot_df[order(plot_df$p_value,decreasing = F),]
  plot_df$term_name=paste(plot_df$source,plot_df$term_name,sep=':')
  plot_df$term_name=factor(plot_df$term_name,levels=rev(plot_df$term_name))
  
  bp=ggplot(plot_df, aes(x = logp, y = term_name)) + 
    geom_point(aes(size = term_size), color = 'blue',alpha = 0.7) +
    labs(x='-log10(p)',y='',title=trait)+
    theme_bw()+
    theme(plot.title=element_text(hjust=0.5),
          axis.text=element_text(size=14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          # legend.position = "bottom"
          )
  #bp
  #ggsave(paste0(save_path,trait,'_bubble.png'),bp,width=11,height=10)#T50
  ggsave(paste0(save_path,trait,'_bubble.png'),bp,width=10.5,height=8)
  
}


