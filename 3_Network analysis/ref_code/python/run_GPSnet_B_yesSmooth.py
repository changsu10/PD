import os

cell_f='../../lesion_data/B_cells_filtered.csv'
save_dir='../Raw_module/B_cells_filtered/'

step=10000
start_index_list=[0,10000,20000,30000,50000]
deg_p_list=[0.05,0.01,0.001]
conn_p_list=[0.05,0.01]
is_smoothing='yes'

for deg_p in deg_p_list:
	for conn_p in conn_p_list:
		for start_index in start_index_list:
			cmd=f'python GPSnet_module_generation_mz.py {cell_f} {save_dir} {start_index} {step} {deg_p} {conn_p} {is_smoothing} &'
			os.system(cmd)