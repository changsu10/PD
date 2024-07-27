library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 1 input: 
      path of files,
      path to save results", call.=FALSE)
}

path=args[1]
save_path=args[2]

traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
         'pigd_scores','tremor_scores','moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
         'gds', 'stai',#mood
         'scopa',#Autonomic
         'ess','rem',#sleep
         'gco',#global
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)

all_common_groups=c('proteinMisfold','oxidativeStress','mitochondrial','neuroinflammation',
         'ferroptosis','axonalTtransport','calciumHomeostasis','gutDysbiosis','Others')

count_df <- data.frame(matrix(0, ncol = length(all_common_groups) + 1, nrow = 0))
colnames(count_df) <- c("Trait", all_common_groups)
count_df$Trait <- as.character(count_df$Trait)

for (trait in traits){
  f=paste0(path,trait,'_enrichedPathway_filtered.tsv')
  if (file.exists(f) && file.info(f)$size != 0 && readLines(f, n = 1) != "\"\"") {
    df=read.table(f,header=TRUE,sep='\t')
    temp_count <- df %>%
      count(common_group) %>%
      spread(common_group, n, fill = 0)
    
    temp_count$Trait <- trait
    temp_count$Trait <- as.character(trait)
    
    missing_groups <- setdiff(all_common_groups, names(temp_count))
    temp_count[missing_groups] <- 0
    temp_count <- temp_count %>%
      dplyr::select(Trait, all_of(all_common_groups))
    
    count_df <- bind_rows(count_df, temp_count)
  }
}

write.table(count_df,paste0(save_path,'pathway_group_count.tsv'),sep='\t',row.names = F)
