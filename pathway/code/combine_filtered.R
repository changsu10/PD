### combine GO, REAC, KEGG filtered into one dataframe per trait
### and map theri level to common categories
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 2 input: 
      save path，
      trait list", call.=FALSE)
}

save_path=args[1]
traits <- as.vector(strsplit(args[2], ",")[[1]])

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
  reac_f=paste0('../REAC/',trait,'_enrichedREAC_filtered.tsv')
  kegg_f=paste0('../KEGG/',trait,'_enrichedKEGG_filtered.tsv')

  df=c()
  if (file.exists(go_f)){
    go=read.table(go_f,sep='\t',header=T, quote = "")
    go=go[go$term_size<1000,]
    df=rbind(df,go)
  }
  
  if (file.exists(reac_f)){
    reac=read.table(reac_f,sep='\t',header=T, quote = "")
    df=rbind(df,reac)
  }

  if (file.exists(kegg_f)){
    kegg=read.table(kegg_f,sep='\t',header=T, quote = "")
    df=rbind(df,kegg)
  }

  df=replace_nullLevel_with_termName(df)
  write.table(df,paste0(trait,'_enrichedPathway_filtered.tsv'),sep='\t',row.names = F)
}

