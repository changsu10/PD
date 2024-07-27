library(circlize)
library(RColorBrewer)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 2 input: 
      read path and save path", call.=FALSE)
}

file_path=args[1]#"../result/planB_countfilter100_norm/combined_all_itself/"
save_path=args[2]#"../result/planB_countfilter100_norm/bubble_all_itself/"

traits=c(#'motor',#'updrs2','updrs3','schwab',#motor
         #'pigd_scores','tremor_scores',
         'moca',#'benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
         'gds', 'stai',#mood
         #'scopa',#Autonomic
         'ess','rem',#sleep
         #'gco',#global
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)

trait='abeta_42'
top_num=30
for (trait in traits){
  df=read.table(paste0(file_path,trait,'_enrichedPathway_filtered.tsv'),sep='\t',header=T)
  df=df[,c('term_name','p_value','source','term_size')]
  ####### bubble plot
  plot_df=df
  plot_df=plot_df[order(plot_df$p_value),][1:top_num,]
  plot_df$logp=-log10(plot_df$p_value)
  # write.csv(plot_df$term_name,paste0(save_path,trait,'_T',top_num,'pathway.csv'),row.names = F,col.names=F,quote = F)
  
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


