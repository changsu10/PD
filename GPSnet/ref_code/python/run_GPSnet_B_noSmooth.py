import os

cell_f='../../lesion_data/B_cells_filtered.csv'
save_dir='../Raw_module/B_cells_filtered/'

step=20000
start_index_list=[0,20000,40000,60000,80000]
deg_p_list=[0.05,0.01,0.001]
conn_p_list=[0.05,0.01]
is_smoothing='no'

for deg_p in deg_p_list:
	for conn_p in conn_p_list:
		for start_index in start_index_list:
			cmd=f'python GPSnet_module_generation_mz.py {cell_f} {save_dir} {start_index} {step} {deg_p} {conn_p} {is_smoothing} &'
			os.system(cmd)