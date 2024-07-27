library('rrvgo')
setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/code/')
GO_result_path='../result/planB_countfilter100_norm/GO/'
save_path='../result_rrvgo/'
traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
         'pigd_scores','tremor_scores','moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
         'gds', 'stai',#mood
         'scopa',#Autonomic
         'ess','rem',#sleep
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)

scatterPlot <- function(simMatrix,reducedTerms,algorithm=c("pca", "umap"),
                        onlyParents=FALSE,addLabel=TRUE,labelSize=3) {
  if(!all(sapply(c("ggplot2", "ggrepel", "umap"), requireNamespace, quietly=TRUE))) {
    stop("Packages ggplot2, ggrepel, umap and/or its dependencies not available. ",
         "Consider installing them before using this function.", call.=FALSE)
  }
  
  if(onlyParents){
    x <- as.data.frame(table(reducedTerms$parentTerm))
    reducedTerms <- reducedTerms[reducedTerms$term == reducedTerms$parentTerm, ]
    simMatrix <- simMatrix[reducedTerms$go, reducedTerms$go]
    reducedTerms[, 'score'] <- x$Freq[match(reducedTerms$term, x$Var1)]
  }
  
  x <- switch(match.arg(algorithm),
              pca =cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)$points,
              umap=umap::umap(as.matrix(as.dist(1-simMatrix)))$layout)
  
  df <- cbind(as.data.frame(x),
              reducedTerms[match(rownames(x), reducedTerms$go), c("term", "parent", "parentTerm", 'score')])
  
  p <-
    ggplot2::ggplot(df, ggplot2::aes(x=V1, y=V2, color=parentTerm,size=score)) +
    ggplot2::geom_point(alpha=.5) +
    ggplot2::scale_color_discrete(guide="none") +
    ggplot2::scale_size_continuous(guide="none", range=c(0, 25)) +
    ggplot2::scale_x_continuous(name="") +
    ggplot2::scale_y_continuous(name="") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank())
  
  if(addLabel) {
    p + ggrepel::geom_label_repel(ggplot2::aes(label=parentTerm),
                                  data=subset(df, parent == rownames(df)),
                                  box.padding=grid::unit(1, "lines"), size=labelSize,
                                  max.overlaps=Inf)
  } else {
    p
  }
}

for (trait in traits){
  go_f=paste0(GO_result_path,trait,'_enrichedGO_filtered.tsv')
  go_df=read.csv(go_f,sep='\t',header=T,col.names=c('term_id','term_name','p_value','source','term_size','ancestor_name'))
  scores <- setNames(-log10(go_df$p_value), go_df$term_id)
  
  ontlogy=c("BP", "MF", "CC")
  for (ont in ontlogy){
    # calculate similarity
    simMatrix <- calculateSimMatrix(go_df$term_id,
                                    orgdb="org.Hs.eg.db",
                                    ont=ont,
                                    method="Rel")
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")
    # plot
    png(filename = paste0(save_path,trait,"_heatmap_",ont,".png"), width = 2200, height = 1600,res=300)
    heatmapPlot(simMatrix,
                reducedTerms,
                annotateParent=TRUE,
                annotationLabel="parentTerm",
                fontsize=6)
    dev.off()
    
    png(filename = paste0(save_path,trait,"_scatter_",ont,".png"), width = 3000, height = 3000,,res=200)
    scatterPlot(simMatrix, reducedTerms)
    dev.off()
    
    png(filename = paste0(save_path,trait,"_treemap_",ont,".png"), width = 2000, height = 2000,,res=200)
    treemapPlot(reducedTerms)
    dev.off()
    
    png(filename = paste0(save_path,trait,"_wordcloud_",ont,".png"), width = 800, height = 800,,res=200)
    wordcloudPlot(reducedTerms, min.freq=1, colors="black")
    dev.off()
    
  }
  print(paste('finish',trait))
}
