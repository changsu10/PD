### KEGG
library(KEGGREST)
library(gprofiler2)
library(data.table)
#library(KEGGgraph)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Need 4 input: 
      GPSnet_result path, 
      save path,
      trait list,
      gpsnet_result_suffix", call.=FALSE)
}

gpsnet_result_path=args[1]
save_path=args[2]
traits <- as.vector(strsplit(args[3], ",")[[1]])
gpsnet_result_suffix = args[4]

queryKEGG <- function(query) {
  result <- NULL
  max_attempts <- 5
  attempt <- 1
  while(is.null(result) && (attempt <= max_attempts)) {
    if (attempt >1){
      message(paste("Attempt", attempt, "of", max_attempts))
    }
    # Try to execute the query
    result <- tryCatch({
      keggGet(paste0('ko',query))[[1]]$CLASS
    }, error = function(e) {
      message("Error occurred:", e$message)
      NULL
    })
    # If the result is still NULL, sleep before retrying
    if (is.null(result)) {
      if (attempt < max_attempts) {
        pause_seconds <- 5 * attempt# Exponential back-off
        message(paste("Waiting for", pause_seconds, "seconds before retrying..."))
        Sys.sleep(pause_seconds)
      }
      attempt <- attempt + 1
    }
  }
  
  return(result)
}

##############################
query <- list()
for (trait in traits) {
  if (trait=='hvlt'){# need to combine 4 into 1
    hvlt=c()
    for (t in c('hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall')){
      file_path <- paste0(gpsnet_result_path, t, gpsnet_result_suffix, '.txt')
      if (file.exists(file_path) && file.info(file_path)$size != 0) {
        hvlt=c(hvlt,as.vector(read.csv(file_path, header = FALSE)$V1))
      }
    }
    hvlt=unique(hvlt)
    query[[trait]]=as.vector(hvlt)
  } else{
    file_path <- paste0(gpsnet_result_path, trait, gpsnet_result_suffix, '.txt')
    # Check if file exists and is not empty
    if (file.exists(file_path) && file.info(file_path)$size != 0) {
      query[[trait]] <- as.vector(read.csv(file_path, header = FALSE)$V1)
    } else {
      query[[trait]] <- vector()  # Creates an empty vector
    }
  }
}

not_found=data.table(triat = character(0), term_id=character(0))
not_found_list=list()
num_notfound=1
for (trait in traits){
  one_trait=query[[trait]]
  one_trait_gostres <- gost(query = one_trait,
                            correction_method='fdr',
                            exclude_iea=TRUE,
                            sources = c("KEGG"),
                            organism = "hsapiens",
                            multi_query = FALSE)
  if (is.null(one_trait_gostres)){
    cat('No enriched pathways in',trait)
    next
  }
  result_table=one_trait_gostres$result
  #result_table=result_table[result_table$term_size>5 & result_table$term_size<1000,]#filter pathways
  terms=result_table$term_id
  terms <- gsub("KEGG:", "", terms)#remove prefix
  
  ## get KEGG Orthology
  level_df <- data.table(term_id = character(0), ancestor_name=character(0))
  rows_list <- list()
  for (i in 1:length(terms)){
    t=terms[i]
    l=queryKEGG(t)
    if (is.null(l)){
      ancestor_name='null'
      not_found_list[[num_notfound]]=list(trait=trait,term_id=t)
      num_notfound=num_notfound+1
    } else{
      seg=strsplit(l,'; ')[[1]]
      ancestor_name=paste(seg,collapse = '|')
    }
    rows_list[[i]] <- list(term_id = t, 
                           ancestor_name=ancestor_name)
  }
  level_df <- rbindlist(list(level_df, rbindlist(rows_list)), use.names = TRUE, fill = TRUE)
  ##
  source=sapply(paste0('KEGG:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$source)
  term_name=sapply(paste0('KEGG:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$term_name)
  term_size=sapply(paste0('KEGG:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$term_size)
  p_value=sapply(paste0('KEGG:',level_df$term_id), function(x) result_table[result_table$term_id==x,]$p_value)
  level_df$source=source
  level_df$term_name=term_name
  level_df$term_size=term_size
  level_df$p_value=p_value
  level_df$term_id=paste0('KEGG:',level_df$term_id)
  level_df=level_df[,c('term_id','term_name','p_value','source','term_size','ancestor_name')]
  write.table(level_df,paste0(save_path,trait,'_enrichedKEGG_filtered.tsv'),sep='\t',row.names = F,quote=F)
  print(paste('! finish',trait,'!'))
}

not_found <- rbindlist(list(not_found, rbindlist(not_found_list)), use.names = TRUE, fill = TRUE)
write.table(not_found,paste0(save_path,'notfound_pathway.tsv'),sep='\t',row.names = F,quote=F)
