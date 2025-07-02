### combine GO, REAC, KEGG filtered into one dataframe per trait
### and map theri level to common categories
# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args)==3) {
#   save_path=args[1]
#   traits <- as.vector(strsplit(args[2], ",")[[1]])
#   databases <- tolower(as.vector(strsplit(args[3], ",")[[1]]))
#   min_size = 5
#   max_size = 1000
# 
# } else if (length(args)==5){
#     save_path=args[1]
#     traits <- as.vector(strsplit(args[2], ",")[[1]])
#     databases <- tolower(as.vector(strsplit(args[3], ",")[[1]]))
#     min_size = as.numeric(args[4])
#     max_size = as.numeric(args[5])
# } else{
#   stop("Need 3 or 5 input: 
#       save path，
#       trait list,
#       databases to combine,
#       min pathway size (optional, default 5),
#       max pathway size (optional, default 1000)", call.=FALSE)
# }

genetic_save_path="~/Desktop/amp_pd_cleaned_results/genetic/pathway.0524/combined_go,reac,kegg/"
genomic_save_path="~/Desktop/amp_pd_cleaned_results/genomic/pathway.0524/combined_go,reac,kegg/"
combine_save_path="~/Desktop/amp_pd_cleaned_results/combined_genetic_geonmic_gene_module_name/pathway.0524/combined_go,reac,kegg/"

all_save_path=c(genetic_save_path,genomic_save_path,combine_save_path)
names=c('genetic','genomic','combined')

traits <- c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
            'pigd_scores','tremor_scores',
            'moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
            'gds', 'stai',#mood
            'scopa',#Autonomic
            'ess','rem',#sleep
            #'gco',#global
            'total_tau','p_tau181p','alpha_syn','abeta_42',#biomarker
            'motor','cognition','mood','sleep'
)
databases <- c('go','reac','kegg')
min_size = 5
max_size = 1000


#####
replace_nullLevel_with_termName=function(df){
  index1=which(df$ancestor_name=='')
  index2=which(df$ancestor_name=='null')
  index=unique(c(index1,index2))
  for (i in index){
    df[i,'ancestor_name']=df[i,'term_name']
  }
  return(df)
}

###
for (idx in c(1,2,3)){
  save_path=all_save_path[idx]
  setwd(save_path)
  summary_list <-c()
  for (trait in traits){
    go_f=paste0('../GO/',trait,'_enrichedGO_filtered.tsv')
    go_no_f=paste0('../GO_nogoxplore_norevigo/',trait,'_enrichedGO_filtered.tsv')
    reac_f=paste0('../REAC/',trait,'_enrichedREAC_filtered.tsv')
    kegg_f=paste0('../KEGG/',trait,'_enrichedKEGG_filtered.tsv')
    
    df <- data.frame()
    num_count <- c(go = 0, reac = 0, kegg = 0, all = 0)
    
    if ('go' %in% databases){
      if (file.exists(go_f)){
        go=read.table(go_f,sep='\t',header=T, quote = "")
        go=go[(go$term_size<=max_size) & (go$term_size>=min_size) ,]
        df=rbind(df,go)
        num_count['go'] <- nrow(go)
      }
    }
    
    if ('reac' %in% databases){
      if (file.exists(reac_f)){
        reac=read.table(reac_f,sep='\t',header=T, quote = "")
        reac=reac[(reac$term_size<=max_size) & (reac$term_size>=min_size),]
        df=rbind(df,reac)
        num_count['reac'] <- nrow(reac)
      }
    }
    
    if ('kegg' %in% databases){
      if (file.exists(kegg_f)){
        kegg=read.table(kegg_f,sep='\t',header=T, quote = "")
        kegg=kegg[(kegg$term_size<=max_size) & (kegg$term_size>=min_size),]
        df=rbind(df,kegg)
        num_count['kegg'] <- nrow(kegg)
      }
    }
    df=replace_nullLevel_with_termName(df)
    num_count['all'] <- nrow(df)
    
    summary_list[[trait]] <- num_count
    
  }
  
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- names(summary_list)
  colnames(summary_df) <- c('GO', 'REAC', 'KEGG', 'All')
  
  write.table(summary_df,paste0('../',names[idx],'_summary_count.tsv'),sep='\t',row.names = T)
  
}

