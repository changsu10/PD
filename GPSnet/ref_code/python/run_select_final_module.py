import os

cutoff_list=[0.005,0.01,0.05,0.1,0.5]
# trait_list=['B_cells','Dendritic_cells','Fibroblasts','Keratinocytes','T_cells','Endothelial_cells',
#             'Plasma_cells','Proliferating_cells','Sweat_gland_Myoepithelial_cells']
trait_list=['Keratinocytes_filtered','T_cells_filtered','Fibroblasts_filtered']
deg_p_list=[0.05,0.01,0.001]
conn_p_list=[0.05,0.01]
is_smoothing='no'
steps='0,20000,40000,60000,80000'
summary_df='../summary.tsv'

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

