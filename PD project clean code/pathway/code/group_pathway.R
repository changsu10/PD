#key words for each common group
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
   
  # ### other groups
  # # metabolism
  # metabolism=c('metabol','Enzyme','energy','Glycolysis','Anabolism','Catabolism'),
  # # genetic information processing
  # geneticInfo=c('genetic','gene','dna','rna','transcript','translat','chromosome',
  #          'genomic','mutation','recombination','regulatory','chromatin','histone'),
  # # response to stress/stimulus
  # stimulusResponse=c('stimulus','stimuli','stress'),
  # # disease
  # disease=c('immune','viral','infecti','disease'),
  # # cell growth, cell cycle, cell structure
  # cellProcess=c('proliferation', 'growth', 'cycle', 'mitotic',
  #          'centrosome', 'vacuole','aggresome','clathrin','membrane'),
  # #signaling pathway
  # signalingPath=c('signal','bind','pathway')
)

num_group=length(groups)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 1 input: 
      version parent path,
      path to save grouped files", call.=FALSE)
}

path=args[1]
save_path=args[2]

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

##
traits=c('updrs1','updrs2','updrs3','updrs4','schwab',#motor
         'pigd_scores','tremor_scores','moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
         'gds', 'stai',#mood
         'scopa',#Autonomic
         'ess','rem',#sleep
         'gco',#global
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)
# traits=c('updrs1')
percent_df=c()
for (trait in traits){
  f=paste0(path,trait,'_enrichedPathway_filtered.tsv')# use all pathways
  
  if (file.exists(f) && file.info(f)$size != 0 && readLines(f, n = 1) != "\"\"") {
    print(trait)
    df=read.table(f,header=TRUE,sep='\t')
    is_group <- matrix(0, nrow = dim(df)[1], ncol = num_group)
    colnames(is_group) = names(groups)
    
    for (i in 1:dim(df)[1]){#for each row
      # all_parents=df[i,]$ancestor_name# use all parent nodes
      all_parents=df[i,]$term_name#use itself
      # all_parents=paste(unique(c(df[i,]$term_name,strsplit(df[i,]$ancestor_name,'\\|')[1][1:2])),collapse = '|')#use itself and the closest parent
      
      matches_any_group <- FALSE
      for (group_name in names(groups)) {
        if (any(sapply(groups[[group_name]], function(keyword) grepl(keyword, all_parents, ignore.case = TRUE)))) {
          is_group[i, which(names(groups) == group_name)] <- 1
          matches_any_group <- TRUE
        }
      }
    }
    
    # convert to dataframe
    is_group_df <- as.data.frame(is_group)
    group_vector <- apply(is_group_df, 1, function(row) {
      non_zero_indices <- which(row != 0)
      if (length(non_zero_indices) == 0) {
        return("Others")
      } else if (length(non_zero_indices) == 1) {
        return(colnames(is_group_df)[non_zero_indices])
      } else {
        return(paste(colnames(is_group_df)[non_zero_indices], collapse = "|"))
      }
    })
    
    save_group_df=df[,c('term_id','term_name','p_value','source','term_size')]
    save_group_df$common_group=group_vector
    
    # save results
    write.table(save_group_df,
                file=paste0(save_path,trait,'_enrichedPathway_filtered.tsv'),
                quote=FALSE,sep='\t',row.names = FALSE)
    
    matches_per_group <- colSums(is_group_df!=0)
    percent_df=rbind(percent_df,matches_per_group/dim(is_group_df)[1]*100)
  } else{
    percent_df=rbind(percent_df,rep(NA,num_group))
    next
  }
  print(paste('finish',trait))
}

percent_df=as.data.frame(percent_df)
rownames(percent_df)=traits
write.csv(percent_df,file=paste0(save_path,"pathway_group_log.csv"),quote=FALSE)




