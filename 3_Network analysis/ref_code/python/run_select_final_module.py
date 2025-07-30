import os

cutoff_list=[0.001,0.003,0.005,0.01]
trait_list=['B_cells_filtered','Endothelial_cells_filtered','Fibroblasts_filtered','Keratinocytes_filtered','T_cells_filtered',
'Myeloid_filtered','Myofibroblast_filtered','Plasma_cells_filtered','Sweat_gland_Myoepithelial_cells_filtered']
deg_p_list=[0.05]
conn_p_list=[0.01]
is_smoothing='yes'
steps='0,10000,20000,30000,40000,50000,60000,70000,80000'
summary_df='../9cellTypes_summary.tsv'

for t in trait_list:
    for cutoff in cutoff_list:
        for deg_p in deg_p_list:
            for conn_p in conn_p_list:
                subfolder=f'degP{deg_p}_conn{conn_p}_{is_smoothing}Smooth'
                save_final_module_path=f'../GPSnet_result/{subfolder}/'
                save_final_module_with_score_path=f'../GPSnet_result_keep_score/{subfolder}/'
                cmd='python select_final_module.py %s %s %s %f %s %s %s' % (t,subfolder,steps,cutoff,save_final_module_path,save_final_module_with_score_path,summary_df)
                os.system(cmd)
    print('finish '+t)

