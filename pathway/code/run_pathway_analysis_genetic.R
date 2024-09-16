setwd('/Users/manage/Desktop/Github/PD/pathway/code/')
library(tibble)

# ---- Set Path ----
combine_database='go,reac,kegg'
gpsnet_result_path0='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/GPSnet/GPSnet_result_final/'
gpsnet_result_score_path0='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/GPSnet/GPSnet_result_keep_score_final/'
parent_save_path='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway'
all_traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
             'pigd_scores','tremor_scores',
             'moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
             'gds', 'stai',#mood
             'scopa',#Autonomic
             'ess','rem',#sleep
             'gco',#global
             'total_tau','p_tau181p','alpha_syn','abeta_42',#biomarker
             'motor','cognition','mood','sleep'
)

trait_list=paste(all_traits,collapse=',')

gpsnet_result_path=gpsnet_result_path0
gpsnet_result_score_path=gpsnet_result_score_path0

summary_list <- list()
gpsnet_result_suffix='.txt'
version = ''
# ---- Run pathway enrichment on 3 databases: GO, REAC, KEGG ----
# GO
save_path=paste0(parent_save_path,version,'/GO/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_goxplore_revigo_filterGO.R",gpsnet_result_path,save_path,trait_list,gpsnet_result_suffix)
system(cmd)

# save enriched number
num_enrichedGO <- sapply(all_traits, function(t) {
  file_path <- paste0(save_path, t, '_enrichedGO_filtered.tsv')
  num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
  num - 1
})
summary_list[[1]] <- num_enrichedGO

# GO - no prune
save_path=paste0(parent_save_path,version,'/GO_nogoxplore_norevigo/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_nogoxplore_norevigo_filterGO.R",gpsnet_result_path,save_path,trait_list,gpsnet_result_suffix)
system(cmd)

# save enriched number
num_enrichedGO_no <- sapply(all_traits, function(t) {
  file_path <- paste0(save_path, t, '_enrichedGO_filtered.tsv')
  num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
  num - 1
})
summary_list[[length(summary_list) + 1]] <- num_enrichedGO_no

# REAC
save_path=paste0(parent_save_path,version,'/REAC/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterREAC.R",gpsnet_result_path,save_path,trait_list,gpsnet_result_suffix)
system(cmd)

# save enriched number
num_enrichedREAC <- sapply(all_traits, function(t) {
  file_path <- paste0(save_path, t, '_enrichedREAC_filtered.tsv')
  num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
  num - 1
})
summary_list[[length(summary_list) + 1]] <- num_enrichedREAC

# KEGG
save_path=paste0(parent_save_path,version,'/KEGG/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterKEGG.R",gpsnet_result_path,save_path,trait_list,gpsnet_result_suffix)
system(cmd)

# save enriched number
num_enrichedKEGG <- sapply(all_traits, function(t) {
  file_path <- paste0(save_path, t, '_enrichedKEGG_filtered.tsv')
  num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
  num - 1
})
summary_list[[length(summary_list) + 1]] <- num_enrichedKEGG


#---- Combine output files into one file ----
combined_save_path=paste0(parent_save_path,version,'/combined_',combine_database,'/')
if (!dir.exists(combined_save_path)) {
  dir.create(combined_save_path, recursive = TRUE)
}
cmd <- paste("Rscript combine_filtered.R",combined_save_path,trait_list,combine_database)
system(cmd)

# save enriched number
num_enriched_combo <- sapply(all_traits, function(t) {
  file_path <- paste0(combined_save_path, t, '_enrichedPathway_filtered.tsv')
  num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
  num - 1
})
summary_list[[length(summary_list) + 1]] <- as.numeric(num_enriched_combo)

# ---- Plot Bubble Plot ----
top_num=100
combined_save_path=paste0(parent_save_path,version,'/combined_',combine_database,'/')
bubble_save_path=paste0(parent_save_path,version,'/bubble_',combine_database,'/T',top_num,'/')
if (!dir.exists(bubble_save_path)) {
  dir.create(bubble_save_path, recursive = TRUE)
}
cmd <- paste("Rscript plot_bubble.R",combined_save_path,bubble_save_path,trait_list,top_num,20,12)#h,w
system(cmd)

top_num=50
bubble_save_path=paste0(parent_save_path,version,'/bubble_',combine_database,'/T',top_num,'/')
if (!dir.exists(bubble_save_path)) {
  dir.create(bubble_save_path, recursive = TRUE)
}
cmd <- paste("Rscript plot_bubble.R",combined_save_path,bubble_save_path,trait_list,top_num,20,12)#h,w
system(cmd)

### Save summary df
summary_df <- as.data.frame(t(do.call(rbind, summary_list)))
colnames(summary_df) <- c('GO(noFilter)','GO_no(noFilter)','REAC(noFilter)','KEGG(noFilter)','GO+REAC+KEGG(filter5-1000)')
summary_df <- rownames_to_column(summary_df, var = "Trait")
write.csv(summary_df,paste0(parent_save_path,'/summary_pathway.csv'),row.names = F, quote = F)

# ---- Collect heatmap data for each database----
heatmap_traits='motor,cognition,gds,stai,scopa,ess,rem'

## GO
database='GO'
pathway_result_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/',database,'/')
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/',database,'/')
if (!dir.exists(heatmap_save_path)) {
  dir.create(heatmap_save_path, recursive = TRUE)
}
cmd=paste('Rscript get_heatmap_table.R',pathway_result_path,heatmap_save_path,heatmap_traits,database)
system(cmd)

## REAC
database='REAC'
pathway_result_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/',database,'/')
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/',database,'/')
if (!dir.exists(heatmap_save_path)) {
  dir.create(heatmap_save_path, recursive = TRUE)
}
cmd=paste('Rscript get_heatmap_table.R',pathway_result_path,heatmap_save_path,heatmap_traits,database)
system(cmd)

## KEGG
database='KEGG'
pathway_result_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/',database,'/')
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/',database,'/')
if (!dir.exists(heatmap_save_path)) {
  dir.create(heatmap_save_path, recursive = TRUE)
}
cmd=paste('Rscript get_heatmap_table.R',pathway_result_path,heatmap_save_path,heatmap_traits,database)
system(cmd)

# ----  Manually order pathways ----
# order_pathway_file: 2 columns csv file with colnames term_name, group

# ----  Plot refined heatmap ----
## GO
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/GO/')
save_name='GO_heatmap_refined'
heatmap_df_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/GO/heatmap_pathway_top20_go.csv'
order_pathway_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/GO/refine_order_pathway_v4.csv'
cmd=paste('Rscript plot_refined_heatmap.R',heatmap_df_file,order_pathway_file,heatmap_traits,heatmap_save_path,save_name)
system(cmd)

## REAC
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/REAC/')
save_name='REAC_heatmap_refined'
heatmap_df_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/REAC/heatmap_pathway_top20_reac.csv'
order_pathway_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/REAC/refine_order_pathway_v2.csv'
cmd=paste('Rscript plot_refined_heatmap.R',heatmap_df_file,order_pathway_file,heatmap_traits,heatmap_save_path,save_name)
system(cmd)

## KEGG
heatmap_save_path=paste0('/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/KEGG/')
save_name='KEGG_heatmap_refined'
heatmap_df_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/KEGG/heatmap_pathway_top20_kegg.csv'
order_pathway_file='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/pathway/heatmap/KEGG/refine_order_pathway_v1.csv'
cmd=paste('Rscript plot_refined_heatmap.R',heatmap_df_file,order_pathway_file,heatmap_traits,heatmap_save_path,save_name)
system(cmd)

