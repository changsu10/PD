library(KEGGREST)
library(gprofiler2)
library(data.table)
library(ReactomeContentService4R)
library("rbioapi")
library(GOxploreR)
library(httr)
library(stringi)
library(XML)
library(readr)
library(GO.db)
library('GOfuncR')
library(ggplot2)
library(dplyr)
library("annotate")    
library(visNetwork) 

gpsnet_result_path='../../GPSnet/planB/matlab/per_trait_filter100_norm/GPSnet_result_final/'
save_path='../result/planB_countfilter100_norm/motor/'

traits=c('motor')

query=list()
motor=c()
for (t in c('updrs2','updrs3','schwab')){
  file_path <- paste0(gpsnet_result_path, t, '.txt')
  if (file.exists(file_path) && file.info(file_path)$size != 0) {
    motor=c(motor,as.vector(read.csv(file_path, header = FALSE)$V1))
  }
}
motor=unique(motor)
query[['motor']]=as.vector(motor)

######## KEGG
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
  result_table=result_table[result_table$term_size>5 & result_table$term_size<1000,]#filter pathways
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

## REAC
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

#### GO
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
  result_table=result_table[result_table$term_size>5,]
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
  
  ## get all ancestors
  ancestor_df=data.table(term_id = character(0), ancestor_name=character(0))
  rows_list <- list()
  for (i in 1:length(revigo_GO)){
    g=revigo_GO[i]
    go_df_hierarchy <- get_parent_nodes(g)
    all_ancestors_name <- paste(go_df_hierarchy$parent_name,collapse ='|')
    if ('|' %in% all_ancestors_name){
      print('improper selection of |')
    }
    rows_list[[i]] <- list(term_id = g, 
                           ancestor_name=all_ancestors_name)
  }
  ancestor_df <- rbindlist(list(ancestor_df, rbindlist(rows_list)), use.names = TRUE, fill = TRUE)
  
  ## append gprofier term_id,	source,	term_name	,term_size,	p_value
  source=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$source)
  term_name=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$term_name)
  term_size=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$term_size)
  p_value=sapply(ancestor_df$term_id, function(x) result_table[result_table$term_id==x,]$p_value)
  
  ancestor_df$source=source
  ancestor_df$term_name=term_name
  ancestor_df$term_size=term_size
  ancestor_df$p_value=p_value
  
  ancestor_df=ancestor_df[,c('term_id','term_name','p_value','source','term_size','ancestor_name')]
  write.table(ancestor_df,paste0(save_path,trait,'_enrichedGO_filtered.tsv'),sep='\t',row.names = F,quote=F)
}

##### Combine
replace_nullLevel_with_termName=function(df){
  index1=which(df$ancestor_name=='')
  index2=which(df$ancestor_name=='null')
  index=unique(c(index1,index2))
  for (i in index){
    df[i,'ancestor_name']=df[i,'term_name']
  }
  return(df)
}

for (trait in traits){
  go_f=paste0(save_path,trait,'_enrichedGO_filtered.tsv')
  reac_f=paste0(save_path,trait,'_enrichedREAC_filtered.tsv')
  kegg_f=paste0(save_path,trait,'_enrichedKEGG_filtered.tsv')

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
  write.table(df,paste0(save_path,trait,'_enrichedPathway_filtered.tsv'),sep='\t',row.names = F)
}

#### Bubble Plot
library(circlize)
library(RColorBrewer)
library(ggplot2)

trait='motor'
top_num=30
df=read.table(paste0(save_path,trait,'_enrichedPathway_filtered.tsv'),sep='\t',header=T)
df=df[,c('term_name','p_value','source','term_size')]
####### bubble plot
plot_df=df
plot_df=plot_df[order(plot_df$p_value),][1:top_num,]
plot_df$logp=-log10(plot_df$p_value)

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

ggsave(paste0(save_path,trait,'_bubble.png'),bp,width=10.5,height=8)

