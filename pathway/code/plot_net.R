## This code will plot a network, where nodes are GPSnet outputs and edges are from PPI.
## Node will be colored by logFC/pathway belongings
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)
library(tidygraph)
library(tidyverse)
library(gprofiler2)
library(RColorBrewer)

# ---- Set File Path ----
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Need 4 input: 
      GPSnet_result path, 
      gpsnet_raw_input_path,
      plot_save_path,
      trait list,
      gpsnet_result_suffix", call.=FALSE)
}

gpsnet_result_path=args[1]
gpsnet_raw_input_path=args[2]
plot_save_path=args[3]
traits <- as.vector(strsplit(args[4], ",")[[1]])
gpsnet_result_suffix = args[5]

# gpsnet_result_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/GPSnet/GPSnet_result_keep_score_final/'
# gpsnet_raw_input_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/GPSnet/lesion_data/'
# plot_save_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/pathway/net_plot/'

# traits <- c('B_cells_filtered', 'Dendritic_cells_filtered', 'Fibroblasts_filtered', 'Endothelial_cells_filtered',
#             'Keratinocytes_filtered', 'Plasma_cells_filtered','Proliferating_cells_filtered',
#             'Sweat_gland_Myoepithelial_cells_filtered','T_cells_filtered')# filenames of GPSnet output files
# traits=c('Keratinocytes_filtered')
#--------------
setwd('/Users/manage/Desktop/Github/PD/pathway/code/')
log2FC_var='log2FoldChange'
gene_var='...1'
pval_var='pvalue'
padj_var='padj'

#--------------
is_label_pathway=FALSE
if (is_label_pathway){
  highlight_pathways <- read_csv('./highlight_pathways.csv')#highlight pathways
}

color_by='logfc'

if (!dir.exists(plot_save_path)) {
  dir.create(plot_save_path, recursive = TRUE)
}

# ---- Read Reference Data ----
edges <- read_csv('../../ref/ppi.csv')
map_ref <- read_csv('../../ref/gene_vocab.csv')
map_ref$ncbi_id=as.character(map_ref$ncbi_id)
G <- graph_from_data_frame(d = edges, directed = FALSE)

# ---- Functions: prep_net, plot_net ----
prep_net=function(t,
                  gpsnet_result_path,gpsnet_result_suffix,
                  gpsnet_raw_input_path,log2FC_var,gene_var,pval_var,padj_var,
                  is_label_pathway){
  # Read module genes
  node_score <- read_csv(paste0(gpsnet_result_path, t, gpsnet_result_suffix, '.txt'), col_names = FALSE)# !! Check if file suffix is txt !!
  colnames(node_score) <- c('ncbi_id', 'gene_confidence_score')
  node_score$gene_confidence_score=as.numeric(node_score$gene_confidence_score)
  nodes_list <- as.character(node_score$ncbi_id)  # Convert to character
  
  # Ensure the nodes in nodes_list are present in the graph
  nodes_list <- nodes_list[nodes_list %in% V(G)$name]
  
  # Select the PPI containing module genes
  H <- induced_subgraph(G, nodes_list)
  
  # Get largest component subgraph
  components <- components(H)
  largest_cc <- which.max(components$csize)
  giantC <- induced_subgraph(H, which(components$membership == largest_cc))
  
  # Get the largest component subgraph edges
  df <- get.data.frame(giantC, what = "edges")
  df <- distinct(df)
  
  # Get largest component subgraph nodes
  nodes_list <- V(giantC)$name
  
  # filter nodes
  node_score <- node_score %>% filter(ncbi_id %in% nodes_list)
  
  # Map gene id to gene symbol
  node_score <- node_score %>% 
    left_join(map_ref %>% select(ncbi_id, symbol), by = 'ncbi_id')
  
  # Append logFC and p_vals
  degs <- read_csv(paste0(gpsnet_raw_input_path, t, '.csv'))# !! Check if raw filenames !!
  node_score <- node_score %>% 
    left_join(degs %>% select({{log2FC_var}}, {{gene_var}}, {{pval_var}}, {{padj_var}}) %>% 
                rename(symbol = {{gene_var}}), by = 'symbol')
  
  # Drop rows with NA values
  node_score <- drop_na(node_score)
  
  # Remove self-loop nodes
  filtered_df <- df %>% filter(from != to) 
  
  common_nodes = intersect(unique(node_score$ncbi_id),unique(c(filtered_df$from,filtered_df$to)))
  node_score = node_score %>% filter(ncbi_id %in% common_nodes)
  filtered_df <- filtered_df %>% filter(from %in% common_nodes & to %in% common_nodes)
  common_nodes <- unique(c(filtered_df$from, filtered_df$to))
  # Filter node_score again to ensure it matches the nodes in filtered_df
  node_score <- node_score %>% filter(ncbi_id %in% common_nodes)
  
  if (is_label_pathway){
    # for each node, find if it belongs to highlight pathways
    pathways=highlight_pathways %>% filter(cell_type==t)
    pathway_gene_map <- list()
    for (pathway_id in pathways$term_id){
      # Map the pathway ID to its associated genes
      gconvert_res <- gconvert(query = pathway_id, organism = "hsapiens", target = "ENTREZGENE_ACC")
      pathway_ncbi_ids <- gconvert_res %>% 
        select(target) %>% 
        pull(target)
      # Store the NCBI IDs in the list with pathway_id as the key
      pathway_gene_map[[pathway_id]] <- pathway_ncbi_ids
    }
    
    node_score <- node_score %>% 
      mutate(pathway = NA_character_)
    
    # Check if each gene belongs to any of the pathways
    for (i in seq_along(node_score$ncbi_id)) {
      gene_id <- as.character(node_score$ncbi_id[i])
      gene_pathways <- pathways$term_id[sapply(pathway_gene_map, function(ids) gene_id %in% ids)]
      if (length(gene_pathways) > 0) {
        if (length(gene_pathways) > 1) {
          #node_score$pathway[i] <- 'Multiple Pathway'
          pathway_name=c()
          for (j in gene_pathways){
            pathway_name=c(pathway_name,pathways$term_name[pathways$term_id == j])
          }
          node_score$pathway[i] <- paste(pathway_name,collapse = '|')
        } else {
          node_score$pathway[i] <- pathways$term_name[pathways$term_id == gene_pathways]
        }
      } else{
        node_score$pathway[i] <- 'Other'
      }
    }
  }
  return(list(filtered_df=filtered_df,node_score=node_score))
}

plot_net=function(t,nodes,edges,is_label_pathway,color_by,log2FC_var){
  # Add edge color
  nodes$ncbi_id <- as.character(nodes$ncbi_id)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)

  if (is_label_pathway){
    # Create a lookup table for pathways
    pathway_lookup <- nodes %>% select(ncbi_id, pathway) %>% deframe()
    # Vectorized operation to determine edge colors
    edges <- edges %>%
      mutate(color = ifelse(pathway_lookup[from] == 'Other' | pathway_lookup[to] == 'Other', '#ececec', 'black'))
    
    # Reorder the levels of the pathway factor dynamically
    unique_pathways <- unique(nodes$pathway)
    reordered_levels <- c(setdiff(unique_pathways, c("Multiple Pathway", "Other")), "Multiple Pathway", "Other")
    nodes$pathway <- factor(nodes$pathway, levels = reordered_levels)
  }
    
  # create graph
  g <- graph_from_data_frame(edges, directed = FALSE,nodes)
  
  # create layout
  layout = create_layout(g, layout = 'stress')
  
  if (color_by=='logfc'){
    # plot1: colored by logFC
    p=ggraph(layout) + 
      geom_edge_link(color='lightgrey',width=0.1) + 
      geom_node_point(aes(size = gene_confidence_score, color = !!sym(log2FC_var))) +
      scale_size_continuous(range = c(5, 10)) + 
      scale_color_gradientn(colours = c('#62aec5','white',"#e64072"),
                            limits = c(-5,5),oob = scales::squish)+
      geom_node_text(aes(label=symbol),size=3)+
      
      theme_graph(base_family = "sans")
    ggsave(paste0(plot_save_path,"network_",t,"_logFC.pdf"), plot = p, device = "pdf", width = 15, height = 15)
  } else if (color_by=='pathway'){
    # plot2: colored by pathway
    p2=ggraph(layout) + 
      geom_edge_link(aes(color = I(color)), width = 0.1) + 
      geom_node_point(aes(size = gene_confidence_score, color = pathway)) +
      scale_color_manual(values=color_vector, breaks = reordered_levels)+
      scale_size_continuous(range = c(5, 10)) + 
      geom_node_text(aes(label=symbol),size=3)+
      theme_graph(base_family = "sans")+
      guides(colour = guide_legend(override.aes = list(size=5)))
    
    ggsave(paste0(plot_save_path,"network_",t,"_pathway.pdf"), plot = p2, device = "pdf", width = 15, height = 15)
    
    # plot3: colored by pathway, all edges grey
    p3=ggraph(layout) + 
      geom_edge_link(color='lightgrey',width=0.1) + 
      geom_node_point(aes(size = gene_confidence_score, color = pathway)) +
      scale_color_manual(values=color_vector, breaks = reordered_levels)+
      scale_size_continuous(range = c(5, 10)) + 
      geom_node_text(aes(label=symbol),size=3)+
      theme_graph(base_family = "sans")+
      guides(colour = guide_legend(override.aes = list(size=5)))
    
    ggsave(paste0(plot_save_path,"network_",t,"_pathway_grey.pdf"), plot = p3, device = "pdf", width = 15, height = 15)
  }
}

for (t in traits) {
  # ---- Prepare Node and Edge ----
  dfs_list=prep_net(t,
                    gpsnet_result_path,gpsnet_result_suffix,
                    gpsnet_raw_input_path,log2FC_var,gene_var,pval_var,padj_var,
                    is_label_pathway)

  # ---- Plot Network ----
  nodes=as.data.frame(dfs_list$node_score)
  edges=dfs_list$filtered_df
  plot_net(t,nodes,edges,is_label_pathway,color_by,log2FC_var)
}

