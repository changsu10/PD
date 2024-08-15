setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/code/')

library('ViSEAGO')
library(data.table)
library(scales)
library(plotly)
library(dplyr)
library(networkD3)
library(stringr)
######## Set variables ########
gene_module_path_genomic='../../GPSnet/planB/matlab/per_trait_filter100_norm/GPSnet_result_keep_score_final/'
gene_module_path_genetic='../../../GPSnet/result_keep_score_final/'

save_path='../result_viseago/'
combined_traits=list('motor'=c('updrs2','updrs3','updrs4','schwab'),
                     'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                     'mood'=c('gds', 'stai'),
                     'autonomic'=c('scopa'),
                     'sleep'=c('ess','rem'),
                     'total_tau'='total_tau',
                     'p_tau181p'='p_tau181p',
                     'alpha_syn'='alpha_syn',
                     'abeta_42'='abeta_42')

######## Generate EntrezGene reference ########
EntrezGene2GO <- function(temp) {# Modify the source Custom2GO function
  gene2go=fread(temp)
  colnames(gene2go) <- c("taxid", "gene_id", "GOID", "evidence", 'Qualifier','GO_term','PubMed',"category")
  
  gene2go[, `:=`(
    taxid = as.character(gene2go$taxid),
    gene_id = as.character(gene2go$gene_id)
  )]
  
  new(
    "genomic_ressource",
    db = "EntrezGene",
    stamp = temp,
    data = gene2go,
    organisms = data.table(taxid=unique(gene2go$taxid))
  )
}

# Use the function with the local file
local_file <- "~/Downloads/gene2go"  # Path to  downloaded file
EntrezGene <- EntrezGene2GO(local_file)

myGENE2GO<-ViSEAGO::annotate(
  "9606",#Taxonomy ID for human
  EntrezGene
)

####### Loop for each trait ########
result_names <- c()
for (comnined_trait in names(combined_traits)){
  # traits in a combined_trait
  traits=unlist(combined_traits[comnined_trait], use.names = FALSE)
  for (path_type in c('genomic','genetic')){
    if (path_type=='genomic'){
      gene_module_path=gene_module_path_genomic
    } else if (path_type=='genetic'){
      gene_module_path=gene_module_path_genetic
    } else{
      print('Wrong path type!!')
    }
    
    files <- paste0(gene_module_path, traits, '.txt')
    # read into one table
    table_list <- lapply(files, function(file) {
      df <- fread(file, sep = ",")
      df
    })
    table <- rbindlist(table_list, use.names = TRUE, fill = TRUE)
    
    # rescale module score to 0-1; the higher score, the more important
    table$V2=1-rescale(table$V2)
    # if has duplicated genes (V1), keep the row with the largest module score (V2)
    table <- table[order(V1, -V2)]
    table <- table[!duplicated(table$V1)]
    data.table::setorder(table,V2)
    
    # perform fgseaMultilevel tests
    BP<-ViSEAGO::runfgsea(
      geneSel=table,
      ont="BP",
      gene2GO=myGENE2GO, 
      method ="fgseaMultilevel",
      params = list(
        scoreType = "pos",
        minSize=5
      )
    )
    result_name <- paste0(comnined_trait,'.',path_type)
    assign(result_name, BP, envir = .GlobalEnv)
    result_names <- c(result_names, result_name)
  }

}  

############## Merge enrichment terms
input_list <- setNames(as.list(result_names), result_names)
BP_sResults <- ViSEAGO::merge_enrich_terms(
  cutoff = 0.1,
  Input = input_list
)

# initialyse 
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

myGOs<-ViSEAGO::compute_SS_distances(myGOs,distance="Wang")
mat<-ifelse(mat < 2, 0, 1)
########### Visualization ###########

# GOterms heatmap, default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(myGOs)
p1=ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)
save_image(p1, file = paste0(save_path,"2all_cluster_heatmap.png"),width=1000,height=2000,scale=10)

# Save GOterms txt
ViSEAGO::show_table(
  Wang_clusters_wardD2,
  paste0(save_path,"2all_cluster_heatmap.txt")
)

# GOterms MDSplot
p2=ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms"
)
save_image(p2, file = paste0(save_path,"2all_mdsplot.png"),width=1200,height=1200,scale=10)

# GOclusters heatmap
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(Wang_clusters_wardD2,distance=c("max", "avg","rcmax", "BMA"))
Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters_wardD2,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)
p3=ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters"
)
save_image(p3, file = paste0(save_path,"2all_cluster_heatmap_group.png"),width=800,height=800,scale=10)

# # Save GOclusters txt
# ViSEAGO::show_table(
#   Wang_clusters_wardD2,
#   paste0(save_path,"2all_GOcluster_heatmap.txt")
# )
############ Sankey
## prepare sankey data
data <- read.table(paste0(save_path,"all_cluster_heatmap.txt"),header=T)

## use which col to calculate percentage
use_term_col='..log10_pvalue'#..log10_pvalue, .NES
use_term_save_name='-log10p'

##
val_columns <- paste0(names(combined_traits),use_term_col)

data_long <- data %>%
  select(GO.cluster, val_columns) %>%
  pivot_longer(cols = -GO.cluster, names_to = "trait", values_to = "value")

data_long <- data_long %>%
  mutate(trait = str_replace(trait, paste0("\\", use_term_col, "$"), ""))

sum_val_table <- data_long %>%
  group_by(GO.cluster, trait) %>%
  summarise(sum_value = sum(`value`, na.rm = TRUE)) %>%
  ungroup()

percentage_table <- sum_val_table %>%
  group_by(trait) %>%
  mutate(percentage = sum_value / sum(sum_value) * 100) %>%
  ungroup()

### Manually extract cluster names
cluster_name <- data.frame(
  GO.cluster = c(2,1,3,4,5,11,9,6,7,8,10,14,15,16,12,13,20,17,18,19,22,24,21,23,
                 25,26,28,27,29,30,31,32,33,45,43,44,39,40,35,36,37,34,38,46,41,42),
  cluster_text = c(
    "protein catabolic process (cl2)",
    "protein metabolic process (cl1) ",
    "regulation of biological process (cl3)",
    "regulation of metabolic process (cl4)" ,
    "gene expression (cl5)",
    "nucleobase-containing compound metabolic process (cl11)" ,
    "cellular macromolecule metabolic process (cl9)",
    "DNA-templated transcription (cl6)" ,
    "regulation of gene expression (cl7)" ,
    "DNA metabolic process (cl8)",
    "DNA repair (cl10)",
    "Wnt signaling pathway (cl14)",
    "signal transduction (cl15)",
    "intracellular signal transduction (cl16)",
    "positive regulation of intracellular signal transd... (cl12)",
    "regulation of signal transduction (cl13)" ,
    "response to organic substance (cl20)",
    "cellular response to cytokine stimulus (cl17)" ,
    "cellular response to endogenous stimulus (cl18)",
    "response to chemical (cl19)",
    "process (cl22)",
    "response to abiotic stimulus (cl24)",
    "response to stimulus (cl21)",
    "response to stimulus (cl23)",
    "animal organ development (cl25)",
    "process (cl26)",
    "neuron projection development (cl28)",
    "cell differentiation (cl27)",
    "regulation of developmental process (cl29)",
    "autophagy (cl30)",
    "cellular component organization (cl31)" ,
    "cellular component organization (cl32)",
    "cellular component assembly (cl33)" ,
    "regulation of molecular function (cl45)" ,
    "cellular metal ion homeostasis (cl43)" ,
    "regulation of biological quality (cl44)" ,
    "protein localization (cl39)",
    "nitrogen compound transport (cl40)" ,
    "mitotic cell cycle phase transition (cl35)",
    "cell migration (cl36)",
    "apoptotic process (cl37)",
    "regulation of cell cycle (cl34)",
    "regulation of cell population proliferation (cl38)",
    "process (cl46)",
    "cell adhesion (cl41)",
    "process (cl42)"
  )
)

merged_data=merge(percentage_table,cluster_name, by = "GO.cluster")

## plot sankey
nodes <- data.frame(name = unique(c(merged_data$trait, merged_data$cluster_text)))
nodes$id <- seq_along(nodes$name) - 1

links <- merged_data %>%
  mutate(source = match(trait, nodes$name) - 1,
         target = match(cluster_text, nodes$name) - 1) %>%
  select(source, target, percentage)

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "source", Target = "target",
                   Value = "percentage", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 20,iterations = 0)
saveWidget(p, file=paste0(save_path,paste0('sankey_',use_term_save_name),".html"))
