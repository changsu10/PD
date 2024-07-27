## CODE
* **preprocess_rna_count.ipynb**\
  Preprocess large raw RNAseq data, and extract basic group/age/sex/race.\
  Output:
  1) PP_rnaseq_filtered.csv
  2) meta.csv
  3) data.csv

* **planB_preprocess_ppmi_clinical.ipynb**\
  Extract clinical data for all samples in PPMI cohort (ignore GWAS data availability & RNA time)\
  Output:
  1) planB_clinical.csv
  2) planB_meta.csv
  3) planB_count.csv

## DATA
**PP_rnaseq_filtered.csv** \
raw RNAseq for samples: 
1) PPMI cohort; 
2) have >1 visit and time gap > 0.5yr.

**meta.csv** \
group/age/sex/race for samples: 
1) PPMI cohort;
2) have >1 visit and time gap > 0.5yr;
3) PD group or healthy control group.

**data.csv** \
RNAseq for samples: 
1) PPMI cohort;
2) have >1 visit and time gap > 0.5yr;
3) PD group or healthy control group.

**planB_clinical.csv** \
clinical scores within 3.5 yrs(<) and static values for all PPMI samples that are:
1) PD;
2) have >1 visit and time gap > 0.5yr. No imputation, only impute MoCA BL by SC and impute scores when calculate additional scores.

**planB_meta.csv** \
similar to meta_filtered.csv but with pd_duration time and ledd and for PD only.

**planB_count.csv** \
similar to count_filtered.csv but for PD only.


## Run
1. run `preprocess_rna_count.ipynb`
2. run `planB_preprocess_ppmi_clinical.ipynb`
3. run `step1_run_lmm_clinical.R`
4. run `step2_run_run_lmm_gene.py`
5. run `step3_cal_beta.R`
6. run `step4_save_per_trait.R`
