library(ReactomeContentService4R)
library(gprofiler2)
library(data.table)
library(rbioapi)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Need 2 input: 
      GPSnet_result path, 
      trait list", call.=FALSE)
}

gpsnet_result_path=args[1]
save_path=args[2]
traits <- as.vector(strsplit(args[3], ",")[[1]])

##############################
query <- list()
for (trait in traits) {
  if (trait=='hvlt'){# need to combine 4 into 1
    hvlt=c()
    for (t in c('hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall')){
      file_path <- paste0(gpsnet_result_path, t, '.txt')
      if (file.exists(file_path) && file.info(file_path)$size != 0) {
        hvlt=c(hvlt,as.vector(read.csv(file_path, header = FALSE)$V1))
      }
    }
    hvlt=unique(hvlt)
    query[[trait]]=as.vector(hvlt)
  } else{
    file_path <- paste0(gpsnet_result_path, trait, '.txt')
    # Check if file exists and is not empty
    if (file.exists(file_path) && file.info(file_path)$size != 0) {
      query[[trait]] <- as.vector(read.csv(file_path, header = FALSE)$V1)
    } else {
      query[[trait]] <- vector()  # Creates an empty vector
    }
  }
}

for (trait in traits){
  one_trait=query[[trait]]
  one_trait_gostres <- gost(query = one_trait,
                            correction_method='fdr',
                            exclude_iea=TRUE,
                            sources = c("REAC"),
                            organism = "hsapiens",
                            multi_query = FALSE)
  if (is.null(one_trait_gostres)){
    cat('No enriched pathways in',trait)
    next
  }
  
  result_table=one_trait_gostres$result
  result_table=result_table[result_table$term_size>5 & result_table$term_size<1000,]#filter pathways
  terms=result_table$term_id
  terms <- gsub("REAC:", "", terms)
  
  ## get top level pathway
  level_df <- data.table(term_id = character(0), ancestor_name=character(0))
  rows_list <- list()
  for (i in 1:length(terms)){
    t=terms[i]
    #l=getPathways(t, top.level = FALSE)$displayName
    l_df=as.data.table(rba_reactome_event_ancestors(event_id=t))
    l=l_df$displayName
    rows_list[[i]] <- list(term_id = t, 
                           ancestor_name=paste(l,collapse = '|'))
  }
  level_df <- rbindlist(list(level_df, rbindlist(rows_list)), use.names = TRUE, fill = TRUE)
  ##
  source=sapply(paste0('REAC:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$source)
  term_name=sapply(paste0('REAC:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$term_name)
  term_size=sapply(paste0('REAC:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$term_size)
  p_value=sapply(paste0('REAC:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$p_value)
  level_df$source=source
  level_df$term_name=term_name
  level_df$term_size=term_size
  level_df$p_value=p_value
  level_df=level_df[,c('term_id','term_name','p_value','source','term_size','ancestor_name')]
  write.table(level_df,paste0(save_path,trait,'_enrichedREAC_filtered.tsv'),sep='\t',row.names = F,quote=F)
  print(paste('! finish',trait,'!'))
}

## visualization code
#exportImage("R-HSA-9701898", output = "diagram", format = "jpg",quality = 10)

