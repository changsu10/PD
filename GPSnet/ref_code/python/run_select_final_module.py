import os

cutoff_list=[0.005,0.01,0.05,0.1,0.5]
# trait_list=['B_cells','Dendritic_cells','Fibroblasts','Keratinocytes','T_cells','Endothelial_cells',
#             'Plasma_cells','Proliferating_cells','Sweat_gland_Myoepithelial_cells']
trait_list=['Keratinocytes']

summary_df='../summary.csv'

for t in trait_list:
    for cutoff in cutoff_list:
        cmd='python select_final_module.py %s %f %s %s %s' % (t,cutoff,'../GPSnet_result/','../GPSnet_result_keep_score/',summary_df)
        os.system(cmd)
    print('finish '+t)

