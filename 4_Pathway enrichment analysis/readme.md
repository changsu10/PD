## Pathway Enrichment Analysis
Related scripts:
* gprofier_goxplore_revigo_filterGO.R: find enriched GO terms and remove duplicated by revigo
* gprofier_filterREAC.R: find enriched REAC terms
* gprofier_filterKEGG.R: find enriched KEGG terms
* combine_filtered.R: combine enriched GO, READ, KEGG terms
* plot_bubble.R: plot bubble plot of top 30 enriched pathways
* run_pathway_analysis.R: run pathway enrichment analysis pipeline

1. Open `run_pathway_analysis.R`. Set working directory to `pathway/code`, and modify paths under `Set Path` section.
2. run `run_pathway_analysis.R`. It will use `gprofier` package to find enriched pathways and plot bubble plots.

* For HS, refer to `run_pathway_analysis_HS.R`

## Group Pathways
### Knowledge Driven -- Pre-defined groups from literature search
Related scripts:
* group_pathway.R: group pathways into predefined groups by keyword searching
* count_groups.R: count the number of groups in each trait
* plot_sankey_diagram.R: plot sankey diagram based on counts

### Data Driven -- Based on GO graph stucture
It's still at exploration stage.
#### ViSEAGO
run_viseago_all.R
#### RRVGO
run_rrvgo_all.R
#### EnrichmentMap in Cytoscape
prepare_enrichmentmap_input.R


## Pathway Visualization
Related scripts:
plot_net.R




