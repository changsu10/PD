setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/code/')

## Run pathway enrichment on gene modules
gpsnet_result_path='../../GPSnet/planB/matlab/per_trait_filter100_norm/GPSnet_result_final/'
parent_save_path='../result/'
version='planB_countfilter100_norm'
suffix='all_itself'
## Run pathway enrichment on 3 databases
# GO
save_path=paste0(parent_save_path,version,'/GO/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_goxplore_revigo_filterGO.R",gpsnet_result_path,save_path)
# cmd <- paste("Rscript gprofier_nogoxplore_norevigo_filterGO.R",gpsnet_result_path,save_path)
system(cmd)
# REAC
save_path=paste0(parent_save_path,version,'/REAC/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterREAC.R",gpsnet_result_path,save_path)
system(cmd)
# KEGG
save_path=paste0(parent_save_path,version,'/KEGG/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterKEGG.R",gpsnet_result_path,save_path)
system(cmd)

## Combine output files into one file
combined_save_path=paste0(parent_save_path,version,'/combined_',suffix,'/')
if (!dir.exists(combined_save_path)) {
  dir.create(combined_save_path, recursive = TRUE)
}
cmd <- paste("Rscript combine_filtered.R",combined_save_path)
system(cmd)
# 
# ## Group pathwasy into common groups
# common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
# cmd <- paste("Rscript group_pathway.R",combined_save_path,common_group_save_path)
# system(cmd)
# 
# ## Generate a count table
# common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
# cmd <- paste("Rscript count_groups.R",common_group_save_path,common_group_save_path)
# system(cmd)

## Plot Sankey diagram
common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
plot_save_path=paste0(parent_save_path,version,'/sankey_',suffix,'/')
if (!dir.exists(plot_save_path)) {
  dir.create(plot_save_path, recursive = TRUE)
}
cmd <- paste("Rscript plot_sankey_diagram.R",common_group_save_path,plot_save_path)
system(cmd)

## Columns normalize heatmap
# heatmap_save_path=paste0(parent_save_path,version,'/heatmap_',suffix,'/')
# if (!dir.exists(heatmap_save_path)) {
#   dir.create(heatmap_save_path, recursive = TRUE)
# }
# cmd <- paste("Rscript plot_heatmap.R",common_group_save_path,heatmap_save_path)
# system(cmd)

## Plot Bubble Plot
combined_save_path=paste0(parent_save_path,version,'/combined_',suffix,'/')
bubble_save_path=paste0(parent_save_path,version,'/bubble_',suffix,'/')
if (!dir.exists(bubble_save_path)) {
  dir.create(bubble_save_path, recursive = TRUE)
}
cmd <- paste("Rscript plot_bubble.R",combined_save_path,bubble_save_path)
system(cmd)



