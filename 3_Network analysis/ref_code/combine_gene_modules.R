## Merge gene modules
library('data.table')
library('dplyr')

combined_traits=list('motor'=c('updrs2','updrs3','schwab'),
                'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                'mood'=c('gds', 'stai'),
                'sleep'=c('ess','rem')
)
gene_module_path='/Users/manage/Desktop/amp_pd_cleaned_results/genetic/GPSnet/GPSnet_result_final/'

# Merge gene modules and save
for (combined_trait in names(combined_traits)) {
  traits <- unlist(combined_traits[combined_trait], use.names = FALSE)
  files <- paste0(gene_module_path, traits, '.txt')
  
  table_list <- lapply(files, function(file) {
    df <- read.csv(file, sep = '\t', header = FALSE)
  })
  
  # Combine the tables for the current combined trait
  combined_table <- rbindlist(table_list, use.names = TRUE, fill = TRUE) %>% distinct()
  write.table(combined_table, paste0(gene_module_path,combined_trait,'.txt'),
            row.names = F, quote = F,col.names = F,
            sep='\t')
}
