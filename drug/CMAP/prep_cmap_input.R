library(dplyr)
library(gprofiler2)
## prepare cmap input up/down data
raw_input_path='~/Desktop/chang_runGPSnet_08-09-2024/DEGs_lesion_vs_normal.psudobulk.default/'
module_file_path='~/Desktop/chang_runGPSnet_08-09-2024/drug_repurposing/degP0.05_conn0.01_yes/'
ref_df=read.csv('~/Desktop/GitHub/PD/ref/gene_vocab.csv')
save_path='~/Desktop/chang_runGPSnet_08-09-2024/drug_repurposing/cmap_code/test_files/'

all_cell_types=c('B_cells','T_cells','Fibroblasts','Keratinocytes')
ref_genes=read.csv('~/Desktop/chang_runGPSnet_08-09-2024/drug_repurposing/cmap_code/all_cmap_gene_ids.csv')
ref_genes=ref_genes$x
#cell_type='B_cells'
for (cell_type in all_cell_types){
  raw_input_df=read.csv(paste0(raw_input_path,cell_type,'_filtered.csv'))
  module_df=read.csv(paste0(module_file_path,cell_type,'_filtered_0.003.txt'))
  
  # convert module ncbi id tp symbol
  module_gene_symbol <- module_df %>%
    left_join(ref_df, by = c("gene" = "ncbi_id"))  %>%
    left_join(raw_input_df, by = c("symbol" = "X"))
  
  up_gene=na.omit(module_gene_symbol[module_gene_symbol$log2FoldChange>0,])$ensembl_id
  down_gene=na.omit(module_gene_symbol[module_gene_symbol$log2FoldChange<0,])$ensembl_id
  
  # convert ensembl to affx id
  up_gene_convert=unique(gconvert(up_gene, organism = "hsapiens", target="AFFY_HG_U133A_2")$target)
  down_gene_convert=unique(gconvert(down_gene, organism = "hsapiens", target="AFFY_HG_U133A_2")$target)
  
  writeLines(up_gene_convert, paste0(save_path,cell_type,"Up.grp"))
  writeLines(down_gene_convert, paste0(save_path,cell_type,"Down.grp"))
  
  #writeLines(intersect(ref_genes,up_gene_convert), paste0(save_path,cell_type,"Up.grp"))
  #print(length(intersect(ref_genes,up_gene_convert)))
  #writeLines(intersect(ref_genes,down_gene_convert), paste0(save_path,cell_type,"Down.grp"))
  #print(length(intersect(ref_genes,down_gene_convert)))
}
