setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/code/')

library('ViSEAGO')
library(data.table)
library(scales)
library(plotly)
##### Set variables
gene_module_path='../../GPSnet/planB/matlab/per_trait_filter100_norm/GPSnet_result_keep_score_final/'
save_path='../result_viseago/'
traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
         ##'pigd_scores','tremor_scores','hvlt',
         'moca','benton','lns','symbol_digit',
         'semantic_fluency',#cognition
         'gds', 'stai',#mood
         'scopa',#Autonomic
         'ess','rem',#sleep
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)

#### Generate EntrezGene reference
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

##### Loop for each trait
for (trait in traits){
  gene_f=paste0(gene_module_path,trait,'.txt')
  table<-data.table::fread(gene_f)
  table$V2=1-rescale(table$V2)#rescale module score to 0-1; the higher score, the more important
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
 
  BP_sResults<-ViSEAGO::merge_enrich_terms(
    cutoff=0.1,
    Input=list(
      condition="BP"
    )
  )
  
  # initialyse 
  myGOs<-ViSEAGO::build_GO_SS(
    gene2GO=myGENE2GO,
    enrich_GO_terms=BP_sResults
  )
  
  myGOs<-ViSEAGO::compute_SS_distances(myGOs,distance="Wang")
  
  ########### Visualization
  # GOterms heatmap, default parameters
  Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(myGOs)
  p=ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOterms"
  )
  save_image(p, file = paste0(save_path,trait,"_cluster_heatmap.png"),scale=10)
  
  # save txt
  ViSEAGO::show_table(
    Wang_clusters_wardD2,
    paste0(save_path,trait,"_cluster_heatmap.txt")
  )
  
  # GOterms MDSplot
  p=ViSEAGO::MDSplot(
    Wang_clusters_wardD2,
    "GOterms"
  )
  save_image(p, file = paste0(save_path,trait,"_mdsplot.png"),scale=10)
  
  
  # GOclusters heatmap
  Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(Wang_clusters_wardD2,distance=c("max", "avg","rcmax", "BMA"))
  Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
    Wang_clusters_wardD2,
    tree=list(
      distance="BMA",
      aggreg.method="ward.D2"
    )
  )
  
  p=ViSEAGO::show_heatmap(
    Wang_clusters_wardD2,
    "GOclusters"
  )
  save_image(p, file = paste0(save_path,trait,"_cluster_heatmap_group.png"),scale=10)
  
  print(paste('Finish',trait,'!!'))

}