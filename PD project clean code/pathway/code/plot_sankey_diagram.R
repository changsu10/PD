# Load necessary libraries
library(dplyr)
library(networkD3)
library(htmlwidgets)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 2 input: 
      read path and save path", call.=FALSE)
}

file_path=args[1]#"../result/planB_countfilter100_norm/group_pathway_all_itself/"
plot_save_path=args[2]#"../result/planB_countfilter100_norm/heatmap_all_itself/"

## read df
df=read.table(paste0(file_path,'pathway_group_count.tsv'),sep='\t',header=TRUE,row.names = 1)
df <- df[, !names(df) %in% "Others"]#remove Others

total_count=c(1782,2741,4411,3436,612,428,1006,238)## GO+REAC+KEGG
total_count=total_count[colSums(df) != 0]
df <- df[, colSums(df) != 0]#remove groups if no trait has pathways in it

## predefine color
node_colors_dict=c(
  'motor'='#ed6d46',
  'updrs2'= '#e15c0c','updrs3'='#ed6d46','schwab'='#f0896a',
  # 'tremor_scores'='#F79C65','pigd_scores'='#f9b186',
  'updrs1'='#B28D72',
  'cognition'='#2e59a7',
  'moca'='#2e59a7','benton'='#527fcf','lns'='#779ad9','hvlt'='#7ea0dc','semantic_fluency'='#bfd0ed',
  'symbol_digit'='#b8caeb',
  
  'mood'='#f9b212',
  'gds'='#fedc5e','stai'='#fac03d',
  
  'autonomic'='#DBC1AA',
  'scopa'='#DBC1AA',
  
  
  'sleep'='#c03b3a',
  'rem'='#cc5756','ess'='#EBB2B5',
  
  'biomarker'='#84c3b7',
  #'gco'='#EB4796','global'='#EB4796',
  'total_tau'='#2a6e3f','p_tau181p'='#3fa65f','alpha_syn'='#7ece96','abeta_42'='#c8ead2',
  
  'axonalTtransport'='#e72382',
  'calciumHomeostasis'='#a6569d',
  'ferroptosis'="#f28cbd",
  'gutDysbiosis'="#b7b2d0",
  'mitochondrial'="#f9cbe1", 
  'neuroinflammation'='#775fa8',
  'oxidativeStress'='#caadd8',
  'proteinMisfold'='#efc0d2'
)

plot_sankey=function(df_normalized,desired_order,node_colors_dict,
                     save_name,plot_save_path){
  # Prepare data for plot
  source <- rep(rownames(df_normalized), times = ncol(df_normalized))#trait
  target <- rep(colnames(df_normalized), each = nrow(df_normalized))#group
  value <- as.vector(as.matrix(df_normalized))
  
  filtered_source <- source[source %in% desired_order]
  filtered_target <- target[source %in% desired_order]
  filtered_value <- value[source %in% desired_order]
  
  nodes <- data.frame(name = unique(c(filtered_source, filtered_target)))
  left_nodes <- nodes %>% filter(name %in% desired_order)
  right_nodes <- nodes %>% filter(!(name %in% desired_order))
  nodes <- bind_rows(left_nodes, right_nodes)
  
  # Create a data frame for links
  links <- data.frame(
    source = match(filtered_source, nodes$name) - 1,
    target = match(filtered_target, nodes$name) - 1,
    value = filtered_value
  )
  
  node_colors <- sapply(nodes$name, function(name) node_colors_dict[name])
  node_colors <- paste0("d3.scaleOrdinal().domain(['", paste(nodes$name, collapse="', '"), "']).range(['", paste(node_colors, collapse="', '"), "'])")
  
  # Create a data frame for links
  links <- data.frame(
    source = match(filtered_source, nodes$name) - 1,
    target = match(filtered_target, nodes$name) - 1,
    value = filtered_value
  )
  
  # set connection colors
  links$group <- as.factor(links$source)
  nodes$group <- as.factor(nodes$name)
  
  domain=paste0('[', paste(sQuote(nodes$name, q = FALSE), collapse = ','), ']')
  colors=as.vector(node_colors_dict[nodes$name])
  range=paste0('[', paste(sQuote(colors, q = FALSE), collapse = ','), ']')
  
  my_color <- paste('d3.scaleOrdinal()', 
                    '.domain(',domain,')', 
                    '.range(',range,')')
  
  #plot
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "source", Target = "target",
                     Value = "value", NodeID = "name", 
                     colourScale=node_colors, LinkGroup="group", NodeGroup="group",
                     sinksRight=FALSE, fontSize = 20,iterations = 0)
  saveWidget(p, file=paste0(plot_save_path,save_name,".html"))
}


##### All traits
desired_order_allTraits <- c('updrs2', 'updrs3', 'schwab', #'pigd_scores','tremor_scores',
                             'updrs1',
                             'moca','benton','lns','hvlt','symbol_digit','semantic_fluency',
                             'gds', 'stai',
                             'scopa',
                             'ess','rem',
                             'total_tau','p_tau181p','alpha_syn','abeta_42')

# colSum norm
df_normalized <- df %>%
  mutate(across(everything(), ~ . / sum(.)))
plot_sankey(df_normalized,desired_order_allTraits,node_colors_dict,
            'allTrait_colSum',plot_save_path)

# total_count norm
df_norm_total <- sweep(df, 2, total_count, "/")
plot_sankey(df_norm_total,desired_order_allTraits,node_colors_dict,
            'allTrait_totalCount',plot_save_path)

## try min-max normalize to 0-1 per row/col
row01_df <- as.data.frame(t(apply(df_norm_total, 1, function(row) {
  min_val <- min(row)
  max_val <- max(row)
  (row - min_val) / (max_val - min_val)
})))
plot_sankey(row01_df,desired_order_allTraits,node_colors_dict,
            'allTrait_totalCount_row0-1',plot_save_path)

### Merged traits
trait_mapping_list=list('motor'=c('updrs4','updrs2','updrs3','schwab'),
                        'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                        'mood'=c('gds', 'stai'),
                        'autonomic'=c('scopa'),
                        'sleep'=c('ess','rem'),
                        'total_tau'=c('total_tau'),
                        'p_tau181p'=c('p_tau181p'),
                        'alpha_syn'=c('alpha_syn'),
                        'abeta_42'=c('abeta_42'))

merge_traits <- function(df, trait_mapping_list) {
  merged_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(merged_df) <- colnames(df)
  for (new_row in names(trait_mapping_list)) {
    merged_df[new_row, ] <- colSums(df[trait_mapping_list[[new_row]],],na.rm=T)
  }
  return(merged_df)
}
desired_order_merged_traits=c('motor','cognition','mood','autonomic','sleep',
                              'total_tau','p_tau181p','alpha_syn','abeta_42')

merged_df <- merge_traits(df, trait_mapping_list)

# colSum norm
merged_df_normalized <- merged_df %>%
  mutate(across(everything(), ~ . / sum(.)))
plot_sankey(merged_df_normalized,desired_order_merged_traits,node_colors_dict,
            'mergedTrait_colSum',plot_save_path)
# total_count norm
merged_df_norm_total <- sweep(merged_df, 2, total_count, "/")
plot_sankey(merged_df_norm_total,desired_order_merged_traits,node_colors_dict,
            'mergedTrait_totalCount',plot_save_path)

## try min-max normalize to 0-1 per row/col
row01_df <- as.data.frame(t(apply(merged_df_norm_total, 1, function(row) {
  min_val <- min(row)
  max_val <- max(row)
  (row - min_val) / (max_val - min_val)
})))
plot_sankey(row01_df,desired_order_merged_traits,node_colors_dict,
            'mregedTrait_totalCount_row0-1',plot_save_path)
