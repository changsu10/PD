library(dbplyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(data.table)
library(rrvgo)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==4){
  GO_result_path=args[1]#paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genomic/pathway/',database,'/')
  save_path=args[2]#paste0('/Users/manage/Desktop/amp_pd/rnaseq/pathway/result_refine/',database,'/')
  all_traits=as.vector(strsplit(args[3], ",")[[1]])
  database=args[4]#'REAC'
  database=toupper(database)
  top_num=20
} else if (length(args)==5){
  GO_result_path=args[1]#paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genomic/pathway/',database,'/')
  save_path=args[2]#paste0('/Users/manage/Desktop/amp_pd/rnaseq/pathway/result_refine/',database,'/')
  all_traits=as.vector(strsplit(args[3], ",")[[1]])
  database=args[4]#'REAC'
  database=toupper(database)
  top_num=as.numeric(args[5])
} else{
  stop("Need 6 input: 
      pathway result path, 
      save path,
      traits,
      database,
      top number", call.=FALSE)
}

# ---- Get initial heatmap data  ----All Traits By Top20 ----
all_top_pathways <- c()
for (one_trait in all_traits) {
  file <- paste0(GO_result_path, one_trait, '_enriched',database,'_filtered.tsv')
  df <- read.csv(file, sep = '\t', header = TRUE, col.names = c('term_id', 'term_name', 'p_value', 'source', 'term_size', 'ancestor_name'))
  df = df[df$source!= 'GO:CC',]# remove GO CC
  df = df[,c('term_id', 'term_name', 'p_value')]
  # Select top pathways
  df = df[order(df$p_value),][1:top_num,]
  all_top_pathways=c(all_top_pathways,df$term_id)
}
all_top_pathways=unique(all_top_pathways)

all_combined_tables <- list()
for (one_trait in all_traits) {
  file <- paste0(GO_result_path, one_trait, '_enriched',database,'_filtered.tsv')
  df <- read.csv(file, sep = '\t', header = TRUE, col.names = c('term_id', 'term_name', 'p_value', 'source', 'term_size', 'ancestor_name'))
  df = df[,c('term_id', 'term_name', 'p_value')]
  
  # Select top pathways
  df = df[df$term_id %in% all_top_pathways, ]
  df$trait=rep(one_trait,dim(df)[1])
  
  # Truncate names to the first 80 characters
  df <- df %>%
    mutate(term_name = ifelse(nchar(term_name) > 85, 
                              paste0(substr(term_name, 1, 80), "...", substr(term_name, nchar(term_name)-4, nchar(term_name))),  
                              term_name))
  # Store the combined table in the list with the name of the combined trait
  all_combined_tables[[one_trait]] <- df
}

go_combined_table <- rbindlist(all_combined_tables, use.names = TRUE, fill = TRUE)
write.csv(go_combined_table,paste0(save_path,'heatmap_pathway_top',top_num,'_',tolower(database),'.csv'),row.names = F)# save heatmap pvalues

# ---- Extract Enriched Trait Count ----
go_combined_table$score=-log10(go_combined_table$p_value)
plot_data0 <- go_combined_table %>%
  dplyr::select(trait, term_name, score) %>%   # Select only relevant columns
  pivot_wider(names_from = term_name, values_from = score,
              values_fn = list(score = mean))   # Reshape to wide format

non_na_counts <- colSums(!is.na(plot_data0))
df_count=data.frame(Pathway=names(non_na_counts[unique(go_combined_table$term_name)]),
                    trait_count=non_na_counts[unique(go_combined_table$term_name)],
                    row.names = NULL)
# write.csv(df_count,
#           paste0(save_path,'pathway_nonNA_trait.csv'),
#           row.names = F,
#           quote = F)
write.table(df_count,paste0(save_path,'pathway_name_count_top',top_num,'_',tolower(database),'.csv'),row.names = F,col.names = F,sep=',')# save unique pathways

# ---- If GO, Use RRVGO to help with clustering ----
if (database=='GO'){
  #go_combined_table=read.csv(paste0(save_path,'heatmap_pathway_top',top_num,'_',tolower(database),'.csv'))
  all_pathway=go_combined_table$term_id
  
  for (ont in c('BP','MF')){
    simMatrix <- calculateSimMatrix(unique(all_pathway),
                                    orgdb="org.Hs.eg.db",
                                    ont=ont,#BP or MF
                                    method="Rel")
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    threshold=0.7,#larger threshold results in fewer clusters
                                    orgdb="org.Hs.eg.db")
    
    # Customize clustering 
    phm=heatmapPlot(simMatrix,
                    reducedTerms,
                    annotateParent=TRUE,
                    annotationLabel="parentTerm",
                    fontsize=6)
    
    row_hclust <- as.hclust(phm$tree_row)
    clusters_row_k <- cutree(row_hclust, k = 10)  # 10 clusters
    customize_clusters = as.data.frame(clusters_row_k)
    colnames(customize_clusters) = c('cluster')
    customize_clusters = rownames_to_column(customize_clusters, var = "go")
    customize_clusters = merge(customize_clusters, reducedTerms[,c('go','score')], by = "go")
    
    # Merge RRVGO clustering results to previous enriched pathway results tables
    colnames(customize_clusters)=c('term_id','cluster','rrvgo_score')
    merged_data <- merge(go_combined_table, customize_clusters, by = "term_id")# merge with customized clustering results
    write.csv(merged_data,paste0(save_path,'rrvgo_cluster_',ont,'.csv'),row.names = F)
  }
}