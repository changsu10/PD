#### pathway enrichment analysis for GO
library(gprofiler2)
library(GOxploreR)
library(httr)
library(stringi)
library(XML)
library(readr)
library(GO.db)
library('GOfuncR')
library(data.table)
library(ggplot2)
library(dplyr)
library(data.table)
library("annotate")    
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
    file_path <- paste0(gpsnet_result_path, trait, gpsnet_result_suffix)
    # Check if file exists and is not empty
    if (file.exists(file_path) && file.info(file_path)$size != 0) {
      query[[trait]] <- as.vector(read.csv(file_path, header = FALSE)$V1)
    } else {
      query[[trait]] <- vector()  # Creates an empty vector
    }

  # if (trait=='hvlt'){# need to combine 4 into 1
  #   hvlt=c()
  #   for (t in c('hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall')){
  #     file_path <- paste0(gpsnet_result_path, t, gpsnet_result_suffix, '.txt')
  #     if (file.exists(file_path) && file.info(file_path)$size != 0) {
  #       hvlt=c(hvlt,as.vector(read.csv(file_path, header = FALSE)$V1))
  #     }
  #   }
  #   hvlt=unique(hvlt)
  #   query[[trait]]=as.vector(hvlt)
  # } else{
  #   file_path <- paste0(gpsnet_result_path, trait, gpsnet_result_suffix, '.txt')
  #   # Check if file exists and is not empty
  #   if (file.exists(file_path) && file.info(file_path)$size != 0) {
  #     query[[trait]] <- as.vector(read.csv(file_path, header = FALSE)$V1)
  #   } else {
  #     query[[trait]] <- vector()  # Creates an empty vector
  #   }
  # }
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
  
  ## no pruning
  save_table=result_table[,c('term_id','term_name','p_value','source','term_size')]
  #save_table$ancestor_name=rep(NA,dim(save_table)[1])
  write.table(save_table,paste0(save_path,trait,'_enrichedGO_filtered.tsv'),sep='\t',row.names = F,quote=F)
  
  # save BP separately
  write.table(save_table[save_table$source=='GO:BP',],paste0(save_path_BP,trait,'_enrichedBP_filtered.tsv'),sep='\t',row.names = F,quote=F)

}






