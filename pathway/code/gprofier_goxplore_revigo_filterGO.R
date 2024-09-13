#### pathway enrichment analysis for GO
library(gprofiler2)
library(GOxploreR)
library(httr)
library(stringi)
library(XML)
library(readr)
library(GO.db)
library(GOfuncR)
library(data.table)
library(ggplot2)
library(dplyr)
library(data.table)
library(annotate)    
library(visNetwork) 

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Need 4 input: 
      GPSnet_result path, 
      save path,
      trait list,
      gpsnet_result_suffix", call.=FALSE)
}

gpsnet_result_path=args[1]
save_path=args[2]

save_path_BP = paste0(save_path,'BP/')
if (!dir.exists(save_path_BP)) {
      dir.create(save_path_BP, recursive = TRUE)
}

traits <- as.vector(strsplit(args[3], ",")[[1]])
gpsnet_result_suffix = args[4]
#################################################
run_revigo <- function(go_file) {
  # Read user data from a file
  fileName <- go_file
  userData <- readChar(fileName, file.info(fileName)$size)
  
  # Submit job to Revigo
  res <- httr::POST(
    url = "http://revigo.irb.hr/Revigo",
    body = list(
      cutoff = "0.5",
      valueType = "PValue",
      speciesTaxon = "0",
      measure = "SIMREL",
      goList = userData
    ),
    encode = "form"
  )
  
  # Check if the response is successful
  if (httr::status_code(res) != 200) {
    print(paste("Failed to get a successful response from Revigo for triat",strsplit(go_file,'_')[[1]][2]))
  }
  
  parsed_html <- XML::htmlParse(res)
  
  # Check if the parsing was successful
  if (is.null(parsed_html)) {
    print(paste("Failed to parse HTML content.",strsplit(go_file,'_')[[1]][2]))
  }
  
  data_table <- data.frame()
  for (i in 1:3) {
    table <- try(XML::readHTMLTable(parsed_html, which = i, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(table, "try-error")) {
      data_table <- rbind(data_table, table)
    }
  }
  
  return(data_table)
}

#################################################
query <- list()
for (trait in traits) {
  if (trait=='hvlt'){# need to combine 4 hvlt into 1
    hvlt=c()
    for (t in c('hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall')){
      file_path <- paste0(gpsnet_result_path, t, gpsnet_result_suffix,'.txt')
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

for (trait in traits){
  one_trait=query[[trait]]
  one_trait_gostres <- gost(query = one_trait,
                            correction_method='fdr',
                            exclude_iea=TRUE,
                            sources = c("GO"),
                            organism = "hsapiens",
                            multi_query = FALSE)
  if (is.null(one_trait_gostres)){
    cat('No enriched pathways in',trait)
    next
  }
  result_table=one_trait_gostres$result
  #result_table=result_table[result_table$term_size>5,]
  goterms=result_table$term_id
  
  ## filter by GOxplore
  a=prioritizedGOTerms(lst=goterms, organism = "Human", sp=TRUE,domain = 'BP')
  b=prioritizedGOTerms(lst=goterms, organism = "Human", sp=TRUE,domain = 'MF')
  c=prioritizedGOTerms(lst=goterms, organism = "Human", sp=TRUE,domain = 'CC')
  filterd_term=c(a$HF,b$HF,c$HF)
  p=c()
  for (f in filterd_term){
    p=c(p,result_table[result_table$term_id==f,]$p_value)
  }
  goxplore_df=data.frame(id=filterd_term,p=p)
  goxplore_file=paste0(save_path,'intermediate_',trait,'_GOxploreR.tsv')
  write.table(goxplore_df,goxplore_file,sep='\t',row.names = F,quote=F,col.names=F)
  
  ### run revigo
  revigo_df=run_revigo(goxplore_file)
  revigo_GO=revigo_df[revigo_df$Representative=='null',]$`Term ID`#select representative GO

  # save if do not need ancestors
  save_df=result_table[result_table$term_id %in% revigo_GO,c('term_id','term_name','p_value','source','term_size')]
  write.table(save_df,paste0(save_path,trait,'_enrichedGO_filtered.tsv'),sep='\t',row.names = F,quote=F)

  # save BP separately
  write.table(save_df[save_df$source=='GO:BP',],paste0(save_path_BP,trait,'_enrichedBP_filtered.tsv'),sep='\t',row.names = F,quote=F)

  # # Comment the following if do not need to get ancestors
  # ## get all ancestors
  # ancestor_df=data.table(term_id = character(0), ancestor_name=character(0))
  # rows_list <- list()
  # for (i in 1:length(revigo_GO)){
  #   g=revigo_GO[i]
  #   go_df_hierarchy <- get_parent_nodes(g)
  #   all_ancestors_name <- paste(go_df_hierarchy$parent_name,collapse ='|')
  #   if ('|' %in% all_ancestors_name){
  #     print('improper selection of |')
  #   }
  #   rows_list[[i]] <- list(term_id = g, 
  #                          ancestor_name=all_ancestors_name)
  # }
  # ancestor_df <- rbindlist(list(ancestor_df, rbindlist(rows_list)), use.names = TRUE, fill = TRUE)
  
  # ## append gprofier term_id,	source,	term_name	,term_size,	p_value
  # source=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$source)
  # term_name=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$term_name)
  # term_size=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$term_size)
  # p_value=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$p_value)
  
  # ancestor_df$source=source
  # ancestor_df$term_name=term_name
  # ancestor_df$term_size=term_size
  # ancestor_df$p_value=p_value
  
  # ancestor_df=ancestor_df[,c('term_id','term_name','p_value','source','term_size','ancestor_name')]
  # write.table(ancestor_df,paste0(save_path,trait,'_enrichedGO_filtered.tsv'),sep='\t',row.names = F,quote=F)
}






