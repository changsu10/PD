## sankey diagram
library(networkD3)
library(tidyr)
library(dplyr)
library(htmlwidgets)
library(tidyverse)

##setwd('/Users/manage/Desktop/amp_pd/pathway/combined_filtered/')
#setwd('/Users/manage/Desktop/amp_pd/rnaseq/pathway/group_pathway/')
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Need 2 input: 
      grouped pathey save path,
       plot save path", call.=FALSE)
}
path=args[1]
setwd(path)

plot_save_path=args[2]#'../sankey_diagram/'

if (!dir.exists(plot_save_path)) {
  dir.create(plot_save_path, recursive = TRUE)
}


traits=c('updrs1','updrs2','updrs3','schwab','updrs4',#motor
         'pigd_scores','tremor_scores','moca','benton','lns','hvlt','symbol_digit','semantic_fluency',#cognition
         'gds', 'stai',#mood
         'scopa',#Autonomic
         'ess','rem',#sleep
         'gco',#global
         'total_tau','p_tau181p','alpha_syn','abeta_42'#biomarker
)

node_colors_dict=c(
  'motor'='#3d9fc0',
  'updrs2'='#e8a5cc','updrs3'='#f282a7','updrs1'='#F1C9E0','updrs4'='#ECB7D6','tremor_scores'='#F79C65','pigd_scores'='#FFD574',
  'cognition'='#84c3b7',
  'moca'='#84c3b7','benton'='#0B9E79','hvlt'='#98ccbb','lns'='#08755A','semantic_fluency'='#B9E0E6',
  'symbol_digit'='#3C967F',
  'mood'='#eaaa60','schwab'='#ad8fd0','scopa'='#caadd8',
  'gds'='#eaaa60','stai'='#9E480B',
  'autonomic'='#DBC1AA',
  'sleep'='#eebbbc',
  'rem'='#eebbbc','ess'='#EBB2B5',
  'biomarker'='#84c3b7',
  'gco'='#EB4796','global'='#EB4796',
  'total_tau'='#84c3b7','abeta_42'='#b8d9c4','p_tau181p'='#d6ead4','alpha_syn'='#eff2db',
  
  'axonalTtransport'='#206491',
  'calciumHomeostasis'='#296cb1',
  'ferroptosis'="#4699c2",
  'gutDysbiosis'="#108b96",
  'mitochondrial'="#55c6d1", 
  'neuroinflammation'='#8db799',
  'oxidativeStress'='#d3e6d3',
  'proteinMisfold'='#3d6036'
  )

############# plot function
plot_sankey=function(df,save_name,node_colors_dict,plot_save_path){
  # df have 2 columns: trait, common_group
  df_cleaned <- df %>%
    group_by(trait, common_group) %>%
    summarise(num = n(), .groups = 'drop')
  df_cleaned=df_cleaned[!is.na(df_cleaned$common_group),]
  
  # link and node
  links=df_cleaned
  colnames(links)=c('source','target','value')
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # set connection colors
  links$group <- as.factor(links$source)
  nodes$group <- as.factor(nodes$name)
  
  domain=paste0('[', paste(sQuote(nodes$name, q = FALSE), collapse = ','), ']')
  colors=as.vector(node_colors_dict[nodes$name])
  range=paste0('[', paste(sQuote(colors, q = FALSE), collapse = ','), ']')
  
  my_color <- paste('d3.scaleOrdinal()', 
                    '.domain(',domain,')', 
                    '.range(',range,')')
  
  #plot
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     colourScale=my_color, LinkGroup="group", NodeGroup="group",
                     sinksRight=FALSE, fontSize = 20)

  # if (!dir.exists(paste0(plot_save_path,save_name))) {
  #   dir.create(paste0(plot_save_path,save_name), recursive = TRUE)
  # }
  saveWidget(p, file=paste0(plot_save_path,save_name,".html"))
}
#############
all_df=c()
for (t in traits){
  f=paste0(t,'_enrichedPathway_filtered.tsv')
  if (file.exists(f) && file.info(f)$size != 0 && readLines(f, n = 1) != "\"\"") {
    df=read.table(f,sep='\t',header=T,quote = "")
    df=df[,c('term_id','p_value','common_group')]
    df$trait=rep(t,dim(df)[1])
    all_df=rbind(all_df,df)
  } else{
    next
  }
}
all_df=all_df[all_df$common_group!='Others',]

# for common_group with |, split
new_df <- all_df %>%
  mutate(common_group = str_split(common_group, "\\|")) %>%
  unnest(common_group)
all_df=new_df

########plot all trait
l1=all_df[,c('trait','common_group')]
plot_sankey(l1,'plot_allTrait',node_colors_dict,plot_save_path)

######## plot all traits, top100 pathways
top100_df <- all_df %>%
  group_by(trait) %>%
  arrange(p_value) %>%
  slice_min(order_by = p_value, n = 100) %>%
  ungroup()

plot_sankey(top100_df,'plot_allTrait_T100',node_colors_dict,plot_save_path)

######## plot all traits, pathways with p values < cutoff
top100_df <- all_df %>%
  group_by(trait) %>%
  filter(p_value < 0.0001) %>%
  arrange(p_value) %>%
  ungroup()

plot_sankey(top100_df,'plot_allTrait_p0.0001',node_colors_dict,plot_save_path)

######## combine traits into motor, cognition, mood, sleep, automatic, mood, biomarker
trait_mapping_list=list('motor'=c('updrs1','updrs4','updrs2','updrs3','schwab','pigd_scores','tremor_scores'),
                        'cognition'=c('moca','benton','lns','hvlt','symbol_digit','semantic_fluency'),
                        'mood'=c('gds', 'stai'),
                        'autonomic'=c('scopa'),
                        'sleep'=c('ess','rem'),
                        'global'=c('gco'),
                        'total_tau'=c('total_tau'),
                        'p_tau181p'=c('p_tau181p'),
                        'alpha_syn'=c('alpha_syn'),
                        'abeta_42'=c('abeta_42')
                        #'biomarker'=c('total_tau','p_tau181p','alpha_syn','abeta_42')
                        )
l2=all_df
mapped_trait=lapply(l2$trait,function(x) for (common_level in names(trait_mapping_list)) {
                                                if (x %in% trait_mapping_list[[common_level]]){return(common_level)}})
l2$trait=unlist(mapped_trait)
plot_sankey(l2,'plot_mergedTrait',node_colors_dict,plot_save_path)

######## plot merged traits, top100 pathways
top100_df <- l2 %>%
  group_by(trait) %>%
  arrange(p_value) %>%
  slice_min(order_by = p_value, n = 100) %>%
  ungroup()

plot_sankey(top100_df,'plot_mergedTrait_T100',node_colors_dict,plot_save_path)


######## plot all traits, pathways with p values < cutoff
top100_df <- l2 %>%
  group_by(trait) %>%
  filter(p_value < 0.0001) %>%
  arrange(p_value) %>%
  ungroup()

plot_sankey(top100_df,'plot_mergedTrait_p0.0001',node_colors_dict,plot_save_path)

