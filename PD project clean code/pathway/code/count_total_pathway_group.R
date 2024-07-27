### get number of pathways in the whole database that containing certain keywords
# library(clusterProfiler)
# library(reactome.db)
# library(org.Hs.eg.db)
# library(KEGGREST)
# library(dplyr)

library(dplyr)
library(stringr)

# Function to read pathway file (assuming the files are in plain text format)
read_pathway_file <- function(file_path) {
  readLines(file_path)
}

go_pathways <- read_pathway_file("/Users/manage/Desktop/amp_pd/rnaseq/pathway/go_name.txt")
reac_pathways <- read_pathway_file("/Users/manage/Desktop/amp_pd/rnaseq/pathway/reac_name.txt")
kegg_pathways <- read_pathway_file("/Users/manage/Desktop/amp_pd/rnaseq/pathway/ko00001.keg")

### 
groups <- list(
  ### PD related groups
  # protein misfolding and PD related protein
  proteinMisfold=c('misfold',#'fold',
                   'chaperone','aggregation',
                   #'stable','stability',
                   'amyloid','conformational','protein homeostasis','proteostasis','prions',
                   'autophagy','ubiquitin','lysosomal',
                   'alpha synuclein','alpha-synuclein','synucleinopathy',
                   'tau','hyperphosphorylation','Lewy'),
  # oxidative stress
  oxidativeStress=c('oxygen','oxidative','oxide','oxida','peroxynitrite','catalase',
                    'peroxidase','redox','free radical','radical','glutathione',
                    'antioxidants','peroxidation','reactive oxygen species'),
  # mitochondrial
  mitochondrial=c('mitochondri','mitophagy','ATP'),
  # neuroinflammation
  neuroinflammation=c('neuroimmun','neuroinflam','neurodegenerat',
                      #'inflam','nerv','immune',
                      'microglia','astrocyte','cytokine','Toll like receptors','tlr',
                      'interleukin','tumor necrosis','chemokine','glial','blood brain'),
  # ferroptosis and cell death
  ferroptosis=c(#'death','kill','aging',
    'ferroptosis','apoptotic',
    'iron','Lipid peroxidation','ferrous'),
  # axonal transport
  axonalTtransport=c('axonal','microtubules','kinesin','dynein','neurofilaments',
                     'anterograde','retrograde','cytoskeletal','neurotrophin','dopaminergic'),
  # calcium homeostasis
  calciumHomeostasis=c('calcium','ryanodine', 'Calmodulin','Inositol trisphosphate'),
  # gut dysbiosis
  gutDysbiosis=c('Microbiota','Gut','Probiotics','Prebiotics','Intestinal','bowel',
                 'dysbiosis','bacteria')#,
)

count_keywords_in_line <- function(line, keywords) {
  sum(sapply(keywords, function(keyword) {
    grepl(keyword, line, ignore.case = TRUE)
  }))
}

### GO
go_results <- data.frame(Line = 1:length(go_pathways))
for (group_name in names(groups)) {
  go_results[[group_name]] <- sapply(go_pathways, function(line) {
    count_keywords_in_line(line, groups[[group_name]])
  })
}
go_count=colSums(go_results)[-1]

### REAC
reac_results <- data.frame(Line = 1:length(reac_pathways))
for (group_name in names(groups)) {
  reac_results[[group_name]] <- sapply(reac_pathways, function(line) {
    count_keywords_in_line(line, groups[[group_name]])
  })
}
reac_count=colSums(reac_results)[-1]

### KEGG
kegg_results <- data.frame(Line = 1:length(kegg_pathways))
for (group_name in names(groups)) {
  kegg_results[[group_name]] <- sapply(kegg_pathways, function(line) {
    count_keywords_in_line(line, groups[[group_name]])
  })
}
kegg_count=colSums(kegg_results)[-1]

count_df=rbind(go_count,reac_count,kegg_count)
count_df=rbind(count_df,colSums(count_df))

write.csv(count_df, "/Users/manage/Desktop/amp_pd/rnaseq/pathway/total_pathway_group_counts.csv", row.names = TRUE)

# > go_count
# proteinMisfold    oxidativeStress      mitochondrial  neuroinflammation        ferroptosis 
# 639                572               2553               1458                377 
# axonalTtransport calciumHomeostasis       gutDysbiosis 
# 150                244                 82 
# > reac_count
# proteinMisfold    oxidativeStress      mitochondrial  neuroinflammation        ferroptosis 
# 31                 52                116                125                 15 
# axonalTtransport calciumHomeostasis       gutDysbiosis 
# 15                  8                  6 
# > kegg_count
# proteinMisfold    oxidativeStress      mitochondrial  neuroinflammation        ferroptosis 
# 1112               2117               1742               1853                220 
# axonalTtransport calciumHomeostasis       gutDysbiosis 
# 263                754                150 

