### Heatmap, column normalize
library(dplyr)
library(ComplexHeatmap)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 2 input: 
      read path and save path", call.=FALSE)
}

file_path=args[1]#"../result/planB_countfilter100_norm/group_pathway_all_itself/"
save_path=args[2]#"../result/planB_countfilter100_norm/heatmap_all_itself/"

df=read.table(paste0(file_path,'pathway_group_count.tsv'),sep='\t',header=TRUE,row.names = 1)
df <- df[, !names(df) %in% "Others"]#remove Others

#proteinMisfold	oxidativeStress	mitochondrial	neuroinflammation	ferroptosis	axonalTtransport	calciumHomeostasis	gutDysbiosis
total_count=c(1782,2741,4411,3436,612,428,1006,238)## GO+REAC+KEGG
# total_count=c(670,	624,	2669,	1583,	392,	165,	252,	88)## GO+REAC
###### All Traits
# divide by colsum
df_norm <- df %>%
  mutate(across(where(is.numeric), ~ . / sum(.)))
png(file=paste0(save_path,'allTrait_divide_colSum.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

# divide by total count
df_norm_total <- sweep(df, 2, total_count, "/")
png(file=paste0(save_path,'allTrait_divide_totalCount.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm_total),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

# z-score
df_norm_z <- df %>%
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
png(file=paste0(save_path,'allTrait_zscore.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm_z),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

## ## try min-max normalize to 0-1 per row/col
row01_df <- as.data.frame(t(apply(df_norm_total, 1, function(row) {
  min_val <- min(row)
  max_val <- max(row)
  (row - min_val) / (max_val - min_val)
})))
png(file=paste0(save_path,'allTrait_totalCount_row0-1.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(row01_df),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

###### Merged Traits
trait_mapping_list=list('motor'=c('updrs1','updrs4','updrs2','updrs3','schwab','pigd_scores','tremor_scores'),
                        'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                        'mood'=c('gds', 'stai'),
                        'autonomic'=c('scopa'),
                        'sleep'=c('ess','rem'),
                        # 'global'=c('gco'),
                        'total_tau'=c('total_tau'),
                        'p_tau181p'=c('p_tau181p'),
                        'alpha_syn'=c('alpha_syn'),
                        'abeta_42'=c('abeta_42')
)
# merge df
merge_traits <- function(df, trait_mapping_list) {
  merged_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(merged_df) <- colnames(df)
  for (new_row in names(trait_mapping_list)) {
    merged_df[new_row, ] <- colSums(df[trait_mapping_list[[new_row]],],na.rm=T)
  }
  return(merged_df)
}
merged_df <- merge_traits(df, trait_mapping_list)

## plot
# divide by colsum
df_norm <- merged_df %>%
  mutate(across(where(is.numeric), ~ . / sum(.)))

png(file=paste0(save_path,'mergedTrait_divide_colSum.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

# divide by total count
df_norm_total <- sweep(merged_df, 2, total_count, "/")
png(file=paste0(save_path,'mergedTrait_divide_totalCount.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm_total),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()


# z-score
df_norm_z <- merged_df %>%
  mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
png(file=paste0(save_path,'mergedTrait_zscore.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(df_norm_z),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()


## try min-max normalize to 0-1 per row/col
row01_df <- as.data.frame(t(apply(df_norm_total, 1, function(row) {
  min_val <- min(row)
  max_val <- max(row)
  (row - min_val) / (max_val - min_val)
})))
png(file=paste0(save_path,'mergedTrait_totalCount_row0-1.png'),width=1000,height=1200,res=200)
h=Heatmap(as.matrix(row01_df),cluster_rows=FALSE,cluster_columns=FALSE)
draw(h)
dev.off()

# col01_df <- df_norm_total %>%
#   mutate(across(everything(), ~ (.-min(.)) / (max(.)-min(.))))
# Heatmap(as.matrix(col01_df),cluster_rows=FALSE,cluster_columns=FALSE)

