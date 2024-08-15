library('rrvgo')
library(networkD3)
library(htmlwidgets)
library(data.table)
library(dplyr)
library(tibble)
library(zeallot)
# ---- Set variables ----
setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/code/')
GO_result_path='../result/planB_countfilter100_norm/GO/'
save_path='../result_rrvgo/'

combined_traits=list('motor'=c('updrs2','updrs3','updrs4','schwab'),
                     'cognition'=c('moca','benton','lns','symbol_digit','semantic_fluency'),
                     'mood'=c('gds', 'stai'),
                     'autonomic'=c('scopa'),
                     'sleep'=c('ess','rem'),
                     'total_tau'='total_tau',
                     'p_tau181p'='p_tau181p',
                     'alpha_syn'='alpha_syn',
                     'abeta_42'='abeta_42')
k=10

save_plot_path=paste0(save_path,'cluster',k,'/')
if (!dir.exists(save_plot_path)) {
  dir.create(save_plot_path, recursive = TRUE)
}

# ---- Combine previous enriched pathway results ----
all_combined_tables <- list()
for (combined_trait in names(combined_traits)) {
  traits <- unlist(combined_traits[combined_trait], use.names = FALSE)
  files <- paste0(GO_result_path, traits, '_enrichedGO_filtered.tsv')
  
  # Read those traits files into one table
  table_list <- lapply(files, function(file) {
    df <- read.csv(file, sep = '\t', header = TRUE, col.names = c('term_id', 'term_name', 'p_value', 'source', 'term_size', 'ancestor_name'))
    df = df[,c('term_id', 'term_name', 'p_value')]
    df
  })
  
  # Combine the tables for the current combined trait
  combined_table <- rbindlist(table_list, use.names = TRUE, fill = TRUE)
  
  # Store the combined table in the list with the name of the combined trait
  all_combined_tables[[combined_trait]] <- combined_table
}

go_combined_table <- rbindlist(all_combined_tables, use.names = TRUE, fill = TRUE)


# ---- Run RRVGO ----
ontlogy='BP'
simMatrix <- calculateSimMatrix(unique(go_combined_table$term_id),
                                orgdb="org.Hs.eg.db",
                                ont=ontlogy,
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.7,#larger threshold results in fewer clusters
                                orgdb="org.Hs.eg.db")

# ---- Extract Dendrogram & Customize clustering ----
phm=heatmapPlot(simMatrix,
                reducedTerms,
                annotateParent=TRUE,
                annotationLabel="parentTerm",
                fontsize=6)
# row_dendrogram <- as.dendrogram(phm$tree_row)
# plot(row_dendrogram, main = "Row Dendrogram")

row_hclust <- as.hclust(phm$tree_row)

# Customize clustering -- Clustering by number of clusters
clusters_row_k <- cutree(row_hclust, k = k)  # k clusters

# # Customize clustering -- Clustering by height
# clusters_row_h <- cutree(row_hclust, h = 6)  # Cut at height h

customize_clusters = as.data.frame(clusters_row_k)
colnames(customize_clusters) = c('cluster')
customize_clusters = rownames_to_column(customize_clusters, var = "go")
customize_clusters = merge(customize_clusters, reducedTerms[,c('go','score')], by = "go")

## Merge RRVGO clustering results to previous enriched pathway results tables
all_go_terms <- bind_rows(all_combined_tables, .id = "trait")
all_go_terms <- all_go_terms %>% rename(go = 'term_id')
merged_data <- merge(all_go_terms, customize_clusters, by = "go")# merge with customized clustering results

# adjust p values
merged_data <- merged_data %>%
  group_by(trait) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

# label clusters by the GO term with highest score. Used in Sankey diagram
cluster_labels <- merged_data %>%
  group_by(cluster) %>%
  slice_max(score, with_ties = FALSE) %>%  # Select the row with the max score, no ties
  select(cluster, term_name) %>%  # Keep only cluster and term_name
  rename(cluster_label = term_name)  # Rename term_name to cluster_label

## save clustering results
save_df = merge(merged_data, cluster_labels, by = "cluster")
write.csv(save_df,file=paste0(save_path,'cluster',k,'_rrvgo_results.csv'),row.names=FALSE)

# ---- Functions: Calculate Normalized Percentage Table ----

#####  Filter pathways options: 1. by pval; 2. by p.adj; 3. by top number
#####  Summarize options: 1. sum p; 2. count N
#####  Normalization options: 1. per trait; 2. per cluster; 3. Jaccard index

## Filter pathways
filter_data <- function(data, 
                        is_filter_by_p = FALSE, 
                        is_filter_by_p_adj = FALSE, 
                        is_filter_by_top = FALSE, 
                        cutoff=1) {
  if (is_filter_by_p) {
    data <- data %>%
      group_by(trait) %>%
      filter(p_value < cutoff)
  } else if (is_filter_by_p_adj) {
    data <- data %>%
      group_by(trait) %>%
      filter(p_adj < cutoff)
  } else if (is_filter_by_top) {
    data <- data %>%
      group_by(trait) %>%
      arrange(p_value) %>%
      slice_head(n = cutoff)
  }
  return(data)
}
## Summarize
sum_data=function(data,
                  is_sum_by_n=FALSE,
                  is_sum_by_p=FALSE,
                  is_sum_by_padj=FALSE){
  if (is_sum_by_n){
    sum_table = data %>%
      distinct(trait, cluster, go) %>%
      group_by(trait, cluster) %>%
      summarise(sum_val = n())
  } else if (is_sum_by_p){
    sum_table = data %>%
      group_by(trait, cluster) %>%
      summarise(sum_val = sum(-log10(p_value)))
  } else if (is_sum_by_padj){
    sum_table = data %>%
      group_by(trait, cluster) %>%
      summarise(sum_val = sum(-log10(p_adj)))
  }
  return(sum_table)
}
#### Normalization
norm_data=function(merged_data,sum_table,
                   is_norm_by_trait=FALSE,
                   is_norm_by_cluster=FALSE,
                   is_norm_by_jaccard=FALSE){
  if (is_norm_by_trait){# how many percentage of pathways in each trait belong to a cluster
    norm_table <- sum_table %>%
      group_by(trait) %>%
      mutate(percentage = sum_val / sum(sum_val) * 100) %>%
      ungroup()
  } else if (is_norm_by_cluster){# how many percentage of pathways in each cluster belong to a trait
    norm_table <- sum_table %>%
      group_by(cluster) %>%
      mutate(percentage = sum_val / sum(sum_val) * 100) %>%
      ungroup()
  } else if (is_norm_by_jaccard){
    trait_cluster_counts <- merged_data %>%
      distinct(trait, cluster, go) %>%
      group_by(trait, cluster) %>%
      summarise(
        count_trait_cluster = n(),
        .groups = 'drop'
      )
    
    trait_counts <- merged_data %>%
      distinct(trait, go) %>%  # Ensure unique pathways per trait
      group_by(trait) %>%
      summarise(
        count_trait = n(),
        .groups = 'drop'
      )
    
    cluster_counts <- merged_data %>%
      distinct(cluster, go) %>%  # Ensure unique pathways per cluster
      group_by(cluster) %>%
      summarise(
        count_cluster = n(),
        .groups = 'drop'
      )
    
    merged_data <- trait_cluster_counts %>%
      left_join(trait_counts, by = "trait") %>%
      left_join(cluster_counts, by = "cluster")
    
    norm_table <- merged_data %>%
      mutate(
        percentage = count_trait_cluster / (count_trait + count_cluster - count_trait_cluster)
      )
    norm_table=norm_table[,c('trait','cluster','percentage')]
  }
  return(norm_table)
}
  
# ---- Function: Plot sankey ----
plot_sankey=function(table,save_path,plot_save_name){
  # prepare data for Sankey diagram
  nodes <- data.frame(name = unique(c(table$trait, table$cluster_label)))
  nodes$id <- seq_along(nodes$name) - 1
  
  links <- table %>%
    mutate(source = match(trait, nodes$name) - 1,
           target = match(cluster_label, nodes$name) - 1) %>%
    select(source, target, percentage)
  # plot
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "source", Target = "target",
                     Value = "percentage", NodeID = "name", 
                     sinksRight=FALSE, fontSize = 20,iterations = 0)
  saveWidget(p, file=paste0(save_path,plot_save_name,".html"))
}

# ---- Run All ----
map <- list(
  filter_by = list(p = c(TRUE, FALSE, FALSE), padj = c(FALSE, TRUE, FALSE), n = c(FALSE, FALSE, TRUE)),
  sum_by = list(p = c(TRUE, FALSE, FALSE), padj = c(FALSE, TRUE, FALSE), n = c(FALSE, FALSE, TRUE)),
  norm_by = list(trait = c(TRUE, FALSE, FALSE), cluster = c(FALSE, TRUE, FALSE), jaccard = c(FALSE, FALSE, TRUE))
)

filter_sum_norm_plot=function(merged_data,filter_by,cutoff,sum_by,norm_by,save_path, plot_save_name, map=map,cluster_labels=cluster_labels){
  # c(is_filter_by_p, is_filter_by_p_adj, is_filter_by_top) %<-% map$filter_by[[filter_by]]
  # c(is_sum_by_p, is_sum_by_padj, is_sum_by_n) %<-% map$sum_by[[sum_by]]
  # c(is_norm_by_trait, is_norm_by_cluster, is_norm_by_jaccard) %<-% map$norm_by[[norm_by]]
  filter_vals <- map$filter_by[[filter_by]]
  is_filter_by_p <- filter_vals[1]
  is_filter_by_p_adj <- filter_vals[2]
  is_filter_by_top <- filter_vals[3]
  
  sum_vals <- map$sum_by[[sum_by]]
  is_sum_by_p <- sum_vals[1]
  is_sum_by_padj <- sum_vals[2]
  is_sum_by_n <- sum_vals[3]
  
  norm_vals <- map$norm_by[[norm_by]]
  is_norm_by_trait <- norm_vals[1]
  is_norm_by_cluster <- norm_vals[2]
  is_norm_by_jaccard <- norm_vals[3]
  
  filtered_data = filter_data(merged_data, is_filter_by_p, is_filter_by_p_adj, is_filter_by_top, cutoff)
  sum_table=sum_data(filtered_data,is_sum_by_n,is_sum_by_p,is_sum_by_padj)
  norm_table=norm_data(merged_data, sum_table,is_norm_by_trait,is_norm_by_cluster,is_norm_by_jaccard)
  
  table=merge(norm_table, cluster_labels, by = "cluster")
  plot_sankey(table,save_path,plot_save_name)
}


filter_by_list=c('p','padj','n')
sum_by_list=c('p','n')#p and p.adj no difference
norm_by_list=c('cluster','jaccard')#trait no pattern

p_cutoff_list=c(0.05,0.01)
padj_cutoff_list=c(0.05,0.01)
n_cutoff_list=c(100)#no difference between different n

for (filter_by in filter_by_list){
  for (sum_by in sum_by_list){
    for (norm_by in norm_by_list){
      if (filter_by=='p'){
        for (cutoff in p_cutoff_list){
          plot_save_name=paste0('sankey_c',k,'_filter.',filter_by,'.',cutoff,'_','sum.',sum_by,'_norm.',norm_by)
          filter_sum_norm_plot(merged_data,filter_by,cutoff,sum_by,norm_by,save_plot_path, plot_save_name, map=map,cluster_labels=cluster_labels)
        }
      } else if (filter_by=='padj'){
        for (cutoff in padj_cutoff_list){
          plot_save_name=paste0('sankey_c',k,'_filter.',filter_by,'.',cutoff,'_','sum.',sum_by,'_norm.',norm_by)
          filter_sum_norm_plot(merged_data,filter_by,cutoff,sum_by,norm_by,save_plot_path, plot_save_name, map=map,cluster_labels=cluster_labels)
        }
      } else if (filter_by=='n'){
        for (cutoff in n_cutoff_list){
          plot_save_name=paste0('sankey_c',k,'_filter.',filter_by,'.',cutoff,'_','sum.',sum_by,'_norm.',norm_by)
          filter_sum_norm_plot(merged_data,filter_by,cutoff,sum_by,norm_by,save_plot_path, plot_save_name, map=map,cluster_labels=cluster_labels)
        }
      }
    }
  }
}

