setwd('/Users/manage/Desktop/Github/PD/pathway/code/')
library(tibble)

# ---- Set Path ----
degP_v=c(0.05,0.01,0.001)
connP_v=c(0.01,0.05)
cutoff_v=c(0.005,0.01,0.05,0.1,0.5)
combine_database_v=c('bp','reac','bp_no','go,reac','go_no,reac')
is_smooth_v=c('noSmooth','yesSmooth')

degP=c(0.05)
connP=c(0.05)
cutoff=c(0.005)
combine_database=c('bp','reac','bp_no','go,reac','go_no,reac')
is_smooth=c('noSmooth')

gpsnet_result_path0='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/GPSnet/GPSnet_result/'
gpsnet_result_score_path0='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/GPSnet/GPSnet_result_keep_score/'
gpsnet_raw_input_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/lesion_data/'
parent_save_path='/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/pathway/KFTB_yenoSmooth/'
trait_list='Fibroblasts_filtered,Keratinocytes_filtered,T_cells_filtered,B_cells_filtered'
trait_list='B_cells_filtered'
all_traits=strsplit(trait_list,',')[[1]]

all_summary_list <- list()
for (degP in degP_v){
  for (connP in connP_v){
    for (is_smooth in is_smooth_v){
      subfolder = paste0('degP',degP,'_conn',connP,'_',is_smooth)
      gpsnet_result_path=paste0(gpsnet_result_path0,subfolder,'/')
      gpsnet_result_score_path=paste0(gpsnet_result_score_path0,subfolder,'/')
      
      for (cutoff in cutoff_v){
        summary_list <- list()
        gpsnet_result_suffix=paste0('_',cutoff)
        version = paste0(subfolder,'_moduleCutoff',cutoff)
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
        num_enrichedBP <- sapply(all_traits, function(t) {
          file_path <- paste0(save_path,'BP/', t, '_enrichedBP_filtered.tsv')
          num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
          num - 1
        })

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
        num_enrichedBP_no <- sapply(all_traits, function(t) {
          file_path <- paste0(save_path,'BP/', t, '_enrichedBP_filtered.tsv')
          num <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1",system(paste("wc -l", file_path), intern = TRUE)))
          num - 1
        })
        
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
        #summary=rbind(summary,num_enrichedREAC)
        
        summary_list[[1]] <- num_enrichedGO
        summary_list[[2]] <- num_enrichedBP
        summary_list[[3]] <- num_enrichedGO_no
        summary_list[[4]] <- num_enrichedBP_no
        summary_list[[5]] <- num_enrichedREAC
        
        for (combine_database in combine_database_v){
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
          #summary=rbind(summary,as.numeric(num_enriched_combo))
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
          
          # ---- Plot Network ----
          net_save_path=paste0(parent_save_path,version,'/net_',combine_database,'/')
          if (!dir.exists(net_save_path)) {
            dir.create(net_save_path, recursive = TRUE)
          }
          cmd <- paste("Rscript plot_net.R",gpsnet_result_score_path,gpsnet_raw_input_path,net_save_path,trait_list,gpsnet_result_suffix)
          system(cmd)
          
        }
        
        summary_df <- as.data.frame(t(do.call(rbind, summary_list)))
        colnames(summary_df) <- c('GO','BP','GO_no','BP_no','REAC',combine_database_v)
        summary_df$degP <- degP
        summary_df$connP <- connP
        summary_df$moduleCut <- cutoff
        summary_df <- rownames_to_column(summary_df, var = "cellType")
        
        all_summary_list[[length(all_summary_list) + 1]] <- summary_df
      }
    }

  }
}

all_summary_df <- do.call(rbind, all_summary_list)
write.csv(all_summary_df,'/Users/manage/Desktop/chang_runGPSnet_08-09-2024/python/pathway/summary.csv')

