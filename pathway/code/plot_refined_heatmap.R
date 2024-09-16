library(dbplyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(data.table)
library(networkD3)
library(htmlwidgets)

### After manually group pathways...
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5){
  stop("Need 5 input: 
      heatmap_df_file, 
      order_pathway_file,
      save_path,
      save_name", call.=FALSE)
}
heatmap_df_file=args[1]#paste0(save_path,'heatmap_pathway_top',top_num,'_',tolower(database),'.csv')
order_pathway_file=args[2]
all_traits=as.vector(strsplit(args[3], ",")[[1]])
save_path=args[4]
save_name=args[5]
# ---- Plot heatmap ----
heatmap_df=read.csv(heatmap_df_file)
heatmap_df$score=-log10(heatmap_df$p_value)

ordered_pathway_df=read.csv(order_pathway_file)## manually change the refined file name
colnames(ordered_pathway_df)=c('term_name','group')
ordered_pathway <- ordered_pathway_df[ordered_pathway_df$term_name != "", 'term_name']

filtered_ordered_pathway <- ordered_pathway[ordered_pathway %in% heatmap_df$term_name]#filter and re-order pathways
heatmap_df=heatmap_df[heatmap_df$term_name %in% filtered_ordered_pathway,]
ordered_pathway=ordered_pathway[ordered_pathway %in% filtered_ordered_pathway]
ordered_pathway_df=ordered_pathway_df[ordered_pathway_df$term_name %in% filtered_ordered_pathway,]

# Prepare heatmap data
plot_data0 <- heatmap_df %>%
  dplyr::select(trait, term_name, score) %>%   # Select only relevant columns
  pivot_wider(names_from = term_name, values_from = score,
              values_fn = list(score = mean))   # Reshape to wide format

plot_data <- as.data.frame(plot_data0)
rownames(plot_data) <- plot_data$trait
plot_data <- plot_data %>% dplyr::select(-trait)
plot_data=plot_data[all_traits,]#re-order rows
plot_data=plot_data[,filtered_ordered_pathway]#re-order columns
plot_matrix <- as.matrix(plot_data)
plot_matrix[is.na(plot_matrix)] <- 0 

# Plot
p=pheatmap(
  mat = t(plot_matrix),
  fontsize = 36,
  fontsize_col = 36,
  cluster_rows = F,  # Do not cluster rows (each row is a file)
  cluster_cols = F,   # Cluster columns based on p_values
  color = colorRampPalette(c("white", '#D53389',"#6D0C68"))(100), # Adjust color palette as needed
  main = "P-value Heatmap"
)
ggsave(paste0(save_path,save_name,'_vertical.pdf'),p,width=28,height=50,dpi=300,limitsize = FALSE)

p=pheatmap(
  mat = plot_matrix,
  fontsize = 36,
  fontsize_col = 36,
  cluster_rows = F,  # Do not cluster rows (each row is a file)
  cluster_cols = F,   # Cluster columns based on p_values
  color = colorRampPalette(c("white", '#D53389',"#6D0C68"))(100), # Adjust color palette as needed
  main = "P-value Heatmap"
)
ggsave(paste0(save_path,save_name,'_horizontal.pdf'),p,width=48,height=30,dpi=300,limitsize = FALSE)

# ---- One group One row: Heatmap ----
groups=unique(ordered_pathway_df$group)
df=merge(x = heatmap_df, y = ordered_pathway_df, by = "term_name") 

# Select the max score (min pval) for each group and each trait
max_score_df <- df %>%
  group_by(group, trait) %>%
  summarize(max_score = max(score, na.rm = TRUE))

group_trait_matrix <- max_score_df %>%
  spread(key = trait, value = max_score)

# Prepare plot data
plot_df=as.data.frame(group_trait_matrix[,-1])
rownames(plot_df)=as.character(group_trait_matrix[,1]$group)
plot_df=plot_df[groups,]#re-order rows
plot_df=plot_df[,all_traits]#re-order cols

plot_matrix=as.matrix(plot_df)
plot_matrix[is.na(plot_matrix)] <- 0 

# Heamtmap
p=pheatmap(plot_matrix,
           fontsize = 36,
           fontsize_col = 36,
           cluster_rows = F, 
           cluster_cols = F, 
           color = colorRampPalette(c("white", '#D53389',"#6D0C68"))(100),
           main = "")
ggsave(paste0(save_path,save_name,'_1group1row.png'),p,width=14,height=20,dpi=300,limitsize = FALSE)
write.csv(plot_matrix,paste0(save_path,save_name,'_1group1row.csv'))

# ---- One group One row: Sankey ----
node_colors_dict=c(
  'motor'='#ed6d46',
  'cognition'='#2e59a7',
  'gds'='#fedc5e','stai'='#fac03d',
  'scopa'='#DBC1AA',
  'rem'='#cc5756','ess'='#EBB2B5')

plot_df_norm=plot_matrix/max(plot_matrix)
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
plot_sankey(plot_df_norm,all_traits,node_colors_dict,save_path,save_name)
