library(circlize)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==4){
  file_path=args[1]
  save_path=args[2]
  traits <- as.vector(strsplit(args[3], ",")[[1]])
  top_num=as.numeric(args[4])
  height=10
  width=12
} else if (length(args)==6){
  file_path=args[1]
  save_path=args[2]
  traits <- as.vector(strsplit(args[3], ",")[[1]])
  top_num=as.numeric(args[4])
  height=as.numeric(args[5])
  width=as.numeric(args[6])
} else{
  stop("Need 6 input: 
      read path, 
      save path,
      traits,
      top number,
      height,
      width", call.=FALSE)
}

for (trait in traits){
  df=read.table(paste0(file_path,trait,'_enrichedPathway_filtered.tsv'),sep='\t',header=T)
  df=df[,c('term_name','p_value','source','term_size')]
  df=df[order(df$p_value),][1:top_num,]
  df$logp=-log10(df$p_value)
  write.csv(df,paste0(save_path,trait,'_bubble_T',top_num,'.csv'),row.names = FALSE)

  df <- df %>%
    mutate(term_name = ifelse(nchar(term_name) > 80, substr(term_name, 1, 80), term_name))
  ####### bubble plot
  plot_df=df
  
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
  ggsave(paste0(save_path,trait,'_bubble_T',top_num,'.png'),bp,width=width,height=height)
  
}


