library(ReactomeContentService4R)
library(gprofiler2)
library(data.table)
library(rbioapi)
library(GOxploreR)
library(httr)
library(stringi)
library(XML)
library(readr)
library(GO.db)
library(GOfuncR)
library(ggplot2)
library(dplyr)
library(annotate)    
library(visNetwork) 

gpsnet_result_path='../../GPSnet/planB/matlab/per_trait_filter100_norm/GPSnet_result_final/'
save_path='../result_for_enrichmentmap/merged_trait/'

# traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
#          'pigd_scores','tremor_scores','moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
#          'gds', 'stai',#mood
#          'scopa',#Autonomic
#          'ess','rem',#sleep
#          'gco',#global
#          'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
# )

trait_mapping_list=list('motor'=c('updrs1','updrs4','updrs2','updrs3','schwab','pigd_scores','tremor_scores'),
                        'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                        'mood'=c('gds', 'stai'),
                        'autonomic'=c('scopa'),
                        'sleep'=c('ess','rem'),
                        'total_tau'=c('total_tau'),
                        'p_tau181p'=c('p_tau181p'),
                        'alpha_syn'=c('alpha_syn'),
                        'abeta_42'=c('abeta_42')
                        )

##############################
query <- list()
for (merged_trait in names(trait_mapping_list)) {
  query_vector=c()
  for (trait in trait_mapping_list[[merged_trait]]) {
    if (trait=='hvlt'){# need to combine 4 into 1
      hvlt=c()
      for (t in c('hvlt_delayed_recall','hvlt_recog_disc_index','hvlt_retention','hvlt_total_recall')){
        file_path <- paste0(gpsnet_result_path, t, '.txt')
        if (file.exists(file_path) && file.info(file_path)$size != 0) {
          hvlt=c(hvlt,as.vector(read.csv(file_path, header = FALSE)$V1))
        }
      }
      hvlt=unique(hvlt)
      # query[[trait]]=as.vector(hvlt)
      query_vector=c(query_vector,hvlt)
    } else{
      file_path <- paste0(gpsnet_result_path, trait, '.txt')
      # Check if file exists and is not empty
      if (file.exists(file_path) && file.info(file_path)$size != 0) {
        query_vector <- c(query_vector,as.vector(read.csv(file_path, header = FALSE)$V1))
      } else {
        # query[[trait]] <- vector()  # Creates an empty vector
      }
    }
  }
  query[[merged_trait]] <- unique(query_vector)
}


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


# for (trait in traits){
for (trait in names(trait_mapping_list)){
  one_trait=query[[trait]]
  #### REAC
  one_trait_gostres_reac <- gost(query = one_trait,
                            correction_method='fdr',
                            exclude_iea=TRUE,
                            evcodes = TRUE,## add this
                            sources = c("REAC"),
                            organism = "hsapiens",
                            multi_query = FALSE)
  if (is.null(one_trait_gostres_reac)){
    cat('No enriched pathways in',trait)
    next
  }
  
  result_table=one_trait_gostres_reac$result
  result_table=result_table[result_table$term_size>5 & result_table$term_size<1000,]#filter pathways
  result_table=result_table[,c("term_id", "term_name", "p_value", "intersection")]
  colnames(result_table) = c("GO.ID", "Description", "p.Val", "Genes")
  result_table$FDR = result_table$p.Val
  result_table$Phenotype = "+1"
  result_table_reac = result_table[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
  # saving the GEM file
  # write.table(result_table, file = paste0(save_path,trait,'_gem.txt'), sep = "\t", quote = F, row.names = F)
  
  print(paste('! finish',trait,'!'))
  
  ## GO
  one_trait_gostres_go <- gost(query = one_trait,
                            correction_method='fdr',
                            exclude_iea=TRUE,
                            evcodes = TRUE,## add this
                            sources = c("GO"),
                            organism = "hsapiens",
                            multi_query = FALSE)
  if (is.null(one_trait_gostres_go)){
    cat('No enriched pathways in',trait)
    next
  }
  
  result_table=one_trait_gostres_go$result
  result_table=result_table[result_table$term_size>5 & result_table$term_size<1000,]
  result_table=result_table[,c("term_id", "term_name", "p_value", "intersection")]
  colnames(result_table) = c("GO.ID", "Description", "p.Val", "Genes")
  result_table$FDR = result_table$p.Val
  result_table$Phenotype = "+1"
  result_table_go = result_table[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
  
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
  
  ##
  result_table_go_filtered=result_table_go[result_table_go$GO.ID %in% revigo_GO, ]
  
  ##
  save_df=rbind(result_table_go_filtered,result_table_reac)
  write.table(save_df, file = paste0(save_path,trait,'_gem.txt'), sep = "\t", quote = F, row.names = F)
}




