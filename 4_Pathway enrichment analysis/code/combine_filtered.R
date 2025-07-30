### combine GO, REAC, KEGG filtered into one dataframe per trait
### and map theri level to common categories
args = commandArgs(trailingOnly=TRUE)

if (length(args)==3) {
  save_path=args[1]
  traits <- as.vector(strsplit(args[2], ",")[[1]])
  databases <- tolower(as.vector(strsplit(args[3], ",")[[1]]))
  min_size = 5
  max_size = 1000

} else if (length(args)==5){
    save_path=args[1]
    traits <- as.vector(strsplit(args[2], ",")[[1]])
    databases <- tolower(as.vector(strsplit(args[3], ",")[[1]]))
    min_size = as.numeric(args[4])
    max_size = as.numeric(args[5])
} else{
  stop("Need 3 or 5 input: 
      save path，
      trait list,
      databases to combine,
      min pathway size (optional, default 5),
      max pathway size (optional, default 1000)", call.=FALSE)
}

setwd(save_path)
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
for (trait in traits){
  go_f=paste0('../GO/',trait,'_enrichedGO_filtered.tsv')
  go_no_f=paste0('../GO_nogoxplore_norevigo/',trait,'_enrichedGO_filtered.tsv')
  reac_f=paste0('../REAC/',trait,'_enrichedREAC_filtered.tsv')
  kegg_f=paste0('../KEGG/',trait,'_enrichedKEGG_filtered.tsv')

  df=c()
  
  if ('go' %in% databases){
    if (file.exists(go_f)){
        go=read.table(go_f,sep='\t',header=T, quote = "")
        go=go[(go$term_size<=max_size) & (go$term_size>=min_size) ,]
        df=rbind(df,go)
    }
  }
  
  if ('go_no' %in% databases){
    if (file.exists(go_no_f)){
      go=read.table(go_no_f,sep='\t',header=T, quote = "")
      go=go[(go$term_size<=max_size) & (go$term_size>=min_size),]
      df=rbind(df,go)
    }
  }
  
  if ('bp' %in% databases){
    bp_f=paste0('../GO/BP/',trait,'_enrichedBP_filtered.tsv')
    if (file.exists(bp_f)){
        go=read.table(bp_f,sep='\t',header=T, quote = "")
        go=go[(go$term_size<=max_size) & (go$term_size>=min_size),]
        df=rbind(df,go)
    }
  }

  if ('bp_no' %in% databases){
    bp_no_f=paste0('../GO_nogoxplore_norevigo/BP/',trait,'_enrichedBP_filtered.tsv')
    if (file.exists(bp_no_f)){
        go=read.table(bp_no_f,sep='\t',header=T, quote = "")
        go=go[(go$term_size<=max_size) & (go$term_size>=min_size),]
        df=rbind(df,go)
    }
  }

  if ('reac' %in% databases){
    if (file.exists(reac_f)){
      reac=read.table(reac_f,sep='\t',header=T, quote = "")
      reac=reac[(reac$term_size<=max_size) & (reac$term_size>=min_size),]
      df=rbind(df,reac)
    }
  }
  
  if ('kegg' %in% databases){
    if (file.exists(kegg_f)){
      kegg=read.table(kegg_f,sep='\t',header=T, quote = "")
      kegg=kegg[(kegg$term_size<=max_size) & (kegg$term_size>=min_size),]
      df=rbind(df,kegg)
    }
 }
  df=replace_nullLevel_with_termName(df)
  write.table(df,paste0(trait,'_enrichedPathway_filtered.tsv'),sep='\t',row.names = F)
}

