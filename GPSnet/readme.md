1. Go the `ref_code`, double check if the `calculate_input` function in `prep_GPSnet_input.py` is the way you want to calculate input (currently it's -log10(p)*beta, modify it if you want to calculate in another way). 
2. Go to `ref_code`, modify `run_prep_GPSnet_input.py` parameters, then run it.
3. Go to the generated folder and the `data` folder inside, manually check the number of input genes in summary_stats.csv. Select the ideal folder (such as pval0.05), copy it to ../ and rename it as `GPSnet_input`.
4. Run the `run_{trait}_xx.m` file. If want to run Matlab on server, can refer the to command line `ref_code/run_all.sh`.
5. Run `run_different_cutoff_keep_score.m`
6. In the generated `GPSnet_result_keep_score` folder, select ideal cutoff results, copy to `../GPSnet_result_keep_score_final` folders
