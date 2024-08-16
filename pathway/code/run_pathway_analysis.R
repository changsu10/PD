setwd('/Users/manage/Desktop/Github/PD/pathway/code/')

# ---- Set Path ----
# gene module path
gpsnet_result_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/GPSnet/DEGs_hs_vs_normal.psudobulk.default/matlab/DEGs_hs_vs_normal.psudobulk.default_norm/GPSnet_result_final/'

# save path
parent_save_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/pathway/'

# save version and suffix
version='DEGs_hs_vs_normal.psudobulk.default'
#suffix='all_itself'

# traits (name of gene module result files)
trait_list='B_cells_filtered,Dendritic_cells_filtered,Fibroblasts_filtered,Keratinocytes_filtered,Plasma_cells_filtered,Proliferating_cells_filtered,Sweat_gland_Myoepithelial_cells_filtered,T_cells_filtered'

# ---- Run pathway enrichment on 3 databases: GO, REAC, KEGG ----
# GO
save_path=paste0(parent_save_path,version,'/GO/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_goxplore_revigo_filterGO.R",gpsnet_result_path,save_path,trait_list)
system(cmd)

# REAC
save_path=paste0(parent_save_path,version,'/REAC/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterREAC.R",gpsnet_result_path,save_path,trait_list)
system(cmd)

# KEGG
save_path=paste0(parent_save_path,version,'/KEGG/')
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
cmd <- paste("Rscript gprofier_filterKEGG.R",gpsnet_result_path,save_path,trait_list)
system(cmd)

#---- Combine output files into one file ----
combine_database='go,reac'
combined_save_path=paste0(parent_save_path,version,'/combined_',combine_database,'/')
if (!dir.exists(combined_save_path)) {
  dir.create(combined_save_path, recursive = TRUE)
}
cmd <- paste("Rscript combine_filtered.R",combined_save_path,trait_list,combine_database)
system(cmd)

# ---- Plot Bubble Plot ----
top_num=100
combined_save_path=paste0(parent_save_path,version,'/combined_',combine_database,'/')
bubble_save_path=paste0(parent_save_path,version,'/bubble_',combine_database,'/T',top_num,'/')
if (!dir.exists(bubble_save_path)) {
  dir.create(bubble_save_path, recursive = TRUE)
}
cmd <- paste("Rscript plot_bubble.R",combined_save_path,bubble_save_path,trait_list,top_num,20,12)#h,w
system(cmd)

# ---- Group pathwasy into common groups ----
# common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
# cmd <- paste("Rscript group_pathway.R",combined_save_path,common_group_save_path)
# system(cmd)
# 
# ## Generate a count table
# common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
# cmd <- paste("Rscript count_groups.R",common_group_save_path,common_group_save_path)
# system(cmd)

# # ---- Plot Sankey diagram ----
# common_group_save_path=paste0(parent_save_path,version,'/group_pathway_',suffix,'/')
# plot_save_path=paste0(parent_save_path,version,'/sankey_',suffix,'/')
# if (!dir.exists(plot_save_path)) {
#   dir.create(plot_save_path, recursive = TRUE)
# }
# cmd <- paste("Rscript plot_sankey_diagram.R",common_group_save_path,plot_save_path)
# system(cmd)





