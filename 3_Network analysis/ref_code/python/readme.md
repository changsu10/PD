## Steps
1. Generate Raw Modules. 
   Assumed folder structure:
   ```
    ├── GPSnet_code/
    └── ref_data/
   ```
   1.1 Copy scripts in `ref_code/python/` into `GPSnet_code`. Copy `ref_data` folder to your path.
   1.2 Run `run_GPSnet_*.py`. Example: `run_GPSnet_9cellTypes.py`
   Output will be:
   ```
   ├── GPSnet_code/
   ├── Raw_Module/
   └── ref_data/
   ```
3. Select final modules.
   Run `run_select_final_module.py`. The parameters in the `run_select_final_module.py` should match the parameters in `run_GPSnet_*.py`.
   Output will be:
   ```
   ├── GPSnet_code/
   ├── GPSnet_result/
   ├── GPSnet_result_keep_score/
   ├── Raw_Module/
   └── ref_data/
   ```

## Notes
1. Difference between `GPSnet_module_generation_mz.py` and `GPSnet_module_generation_mz_limitSize.py` is line120-121. If the raw module size is too big, it will break the loop to save time.
2. Before running, make sure how many combinations you need to run and reserve a server with sufficient cores.
