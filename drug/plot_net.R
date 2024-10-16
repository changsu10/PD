library(igraph)
library(R.matlab)
setwd('~/Desktop/amp_pd/drug_repurposing/')
save_plot_path='./plot/network/'

# set path
drug_df=read.csv('./data/drug_candidates.csv')
drug_df[drug_df$Trait=='non-motor',]$Trait='non_motor'
drug_traits=unique(drug_df$Trait)

traits=list(motor=c('updrs2','updrs3','schwab'),
            cognition=c('moca','benton','lns','symbol_digit','semantic_fluency','hvlt'),
            mood=c('gds', 'stai'),
            sleep=c('ess', 'rem'),
            autonomic= c('scopa'),
            common=c('updrs2','updrs3','schwab','moca','benton','lns','symbol_digit','semantic_fluency','hvlt','gds', 'stai','ess', 'rem','scopa'),
            non_motor=c('moca','benton','lns','symbol_digit','semantic_fluency','hvlt','gds', 'stai','ess', 'rem','scopa'))
module_gene_path='~/Desktop/amp_pd/combined_genetic_geonmic_gene_module/'

all_module_genes=list()
for (i in drug_traits){
  for (j in traits[[i]]){
    all_module_genes[[i]]=c(all_module_genes[[i]],read.csv(paste0(module_gene_path,j,'.txt'),header=F)$V1)
  }
  all_module_genes[[i]]=unique(all_module_genes[[i]])
}


### read ref data
ref_mat<- readMat('./Data_mat/Gene_Drug_10uM.mat')
gene_drug <- ref_mat$Gene.Drug
drug_list <- sapply(ref_mat$Drug.List, function(x) x[[1]][,1])
gene_drug_df <- data.frame(gene_id = gene_drug[, 1], drug = drug_list[gene_drug[,2]])

PPI <- readMat('./Data_mat/Net_PPI.mat')
PPI=PPI$Net
PPI_df=PPI[,c(1,2)]

ref_gene=read.csv('../ref/gene_vocab.csv')

### 
for (i in 1:dim(drug_df)[1]){
  trait=drug_df[i,1]
  drug_name=drug_df[i,2]
  module_genes=all_module_genes[[trait]]
    
  sub_drug_gene_df=gene_drug_df[gene_drug_df$drug==drug_name,]
  colnames(sub_drug_gene_df)=c('ncbi_id','drug')
  
  #connection between module genes and drug genes
  sub_ppi=PPI_df[((PPI_df[,1] %in% module_genes) & (PPI_df[,2] %in% sub_drug_gene_df$ncbi_id)) | 
                   ((PPI_df[,2] %in% module_genes) & (PPI_df[,1] %in% sub_drug_gene_df$ncbi_id)), ]
  sub_ppi=as.data.frame(matrix(sub_ppi,ncol=2))
  
  # Combine the edges
  colnames(sub_ppi)=c('from','to')
  colnames(sub_drug_gene_df)=c('from','to')
  all_edges <- rbind(sub_drug_gene_df, sub_ppi)
  
  # Assign colors based on node type
  nodes <- unique(c(all_edges$from, all_edges$to))
  node_colors <- rep("grey", length(nodes))  # Default color
  
  node_colors[nodes %in% c(sub_ppi[,1],sub_ppi[,2])] <- "#e3eefd"#blue module genes
  node_colors[nodes %in% sub_drug_gene_df$from] <- "#b8d9c4"#green drug targets
  node_colors[nodes %in% sub_drug_gene_df$to] <- "#fbddc1"#drug
  
  ##
  g <- graph_from_data_frame(d = all_edges, directed = FALSE)
  V(g)$color <- node_colors
  new_labels <- V(g)$name
  is_numeric <- grepl("^[0-9]+$", V(g)$name)
  new_labels[is_numeric] <- ref_gene$symbol[match(as.numeric(V(g)$name[is_numeric]), ref_gene$ncbi_id)]
  V(g)$label <- new_labels
  
  # Plot the network
  png(paste0(save_plot_path,trait,"_",drug_name,".png"), width = 10, height = 10, units='in',res=300)
  one_drug=plot(g, vertex.size = 15, 
                vertex.label.cex = 0.7, 
                vertex.color = V(g)$color, 
                edge.color = "grey", 
                main = paste(trait,'-',drug_name))
  legend("topright", legend = c("Module Genes", "Drug Targets", "Drugs"),
         col = c("#e3eefd", "#b8d9c4", "#fbddc1"), pch = 19, cex = 0.8)
  dev.off()
}


