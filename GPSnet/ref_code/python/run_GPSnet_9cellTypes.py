import os

cell_f_list=['./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/B_cells_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Endothelial_cells_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Fibroblasts_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Keratinocytes_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Myeloid_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Myofibroblast_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Plasma_cells_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/Sweat_gland_Myoepithelial_cells_filtered.csv',
'./Outputs_nyu_external_external2_external3.1202/DEGs_lesion_vs_normal.psudobulk/T_cells_filtered.csv']

save_dir_list=['./Raw_module/B_cells_filtered/',
'./Raw_module/Endothelial_cells_filtered/',
'./Raw_module/Fibroblasts_filtered/',
'./Raw_module/Keratinocytes_filtered/',
'./Raw_module/Myeloid_filtered/',
'./Raw_module/Myofibroblast_filtered/',
'./Raw_module/Plasma_cells_filtered/',
'./Raw_module/Sweat_gland_Myoepithelial_cells_filtered/',
'./Raw_module/T_cells_filtered/']

step=10000
start_index_list=[0,10000,20000,30000,40000,50000,60000,70000,80000,90000]
deg_p_list=[0.05]
conn_p_list=[0.01]
is_smoothing='yes'

for i in range(len(cell_f_list)):
	cell_f=cell_f_list[i]
	save_dir=save_dir_list[i]
	for deg_p in deg_p_list:
		for conn_p in conn_p_list:
			for start_index in start_index_list:
				cmd=f'python GPSnet_module_generation_mz_limitSize.py {cell_f} {save_dir} {start_index} {step} {deg_p} {conn_p} {is_smoothing} &'
				os.system(cmd)