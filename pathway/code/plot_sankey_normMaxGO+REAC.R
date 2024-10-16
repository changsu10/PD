library(dbplyr)
library(dplyr)
library(tidyr)
library(data.table)
library(networkD3)
library(htmlwidgets)

go_df_file="~/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/GO/GO_heatmap_refined_1group1row.csv"
reac_df_file="~/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/REAC/REAC_heatmap_refined_1group1row.csv"
save_path="~/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/sankey/"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
all_traits=c('motor','cognition','gds','stai','scopa','rem','ess')
plot_matrix_go=read.csv(go_df_file,row.names = 1)
plot_matrix_reac=read.csv(reac_df_file,row.names = 1)
max_val=max(max(plot_matrix_go),max(plot_matrix_reac))

# ---- Function ----
plot_sankey=function(plot_df,all_traits,node_colors_dict,save_path,save_name){
  # Prepare data for plot
  target <- rep(rownames(plot_df), times = ncol(plot_df))#trait
  source <- rep(colnames(plot_df), each = nrow(plot_df))#group
  value <- as.vector(as.matrix(plot_df))
  
  filtered_source <- source[source %in% all_traits]
  filtered_target <- target[source %in% all_traits]
  filtered_value <- value[source %in% all_traits]
  
  nodes <- data.frame(name = unique(c(filtered_source, filtered_target)))
  left_nodes <- nodes %>% filter(name %in% all_traits)
  right_nodes <- nodes %>% filter(!(name %in% all_traits))
  nodes <- bind_rows(left_nodes, right_nodes)
  
  # Create a data frame for links
  links <- data.frame(
    source = match(filtered_source, nodes$name) - 1,
    target = match(filtered_target, nodes$name) - 1,
    value = filtered_value
  )
  
  node_colors <- sapply(nodes$name, function(name) {
    if (name %in% names(node_colors_dict)) {
      node_colors_dict[name]
    } else {
      "#e5e5e5"
    }
  })
  names(node_colors) <- nodes$name
  
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
                     colourScale=node_colors,LinkGroup="group", NodeGroup="group",
                     sinksRight=FALSE, fontSize = 20,iterations = 0)
  
  saveWidget(p, file=paste0(save_path,save_name,'_sankey',".html"))
}

# ---- Plot ----
node_colors_dict=c('motor'='#ed6d46','cognition'='#2e59a7','gds'='#fedc5e','stai'='#fac03d','scopa'='#DBC1AA','rem'='#cc5756','ess'='#EBB2B5')
plot_sankey(plot_matrix_go/max_val,all_traits,node_colors_dict,save_path,'GO_normMaxGO+REAC')
plot_sankey(plot_matrix_reac/max_val,all_traits,node_colors_dict,save_path,'REAC_normMaxGO+REAC')
