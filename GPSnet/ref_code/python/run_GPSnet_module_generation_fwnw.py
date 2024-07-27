import os
import sys

input_path='../FMNW_interpolation/GPSnet_input/pval0.1/'
DIR_MODULE = '../FMNW_interpolation/Raw_Modules/pval0.1/'
INPUT_PPI_FILE='../ref/ppi.csv'

if not os.path.exists(DIR_MODULE):
    os.makedirs(DIR_MODULE)

file_list=os.listdir(input_path)
file_list=[i for i in file_list if i.endswith('.csv')]

for f in file_list:
    trait=('_').join(f.split('_')[:-2])
    out_path=DIR_MODULE+trait
    # if not os.path.exists(out_path):
    #     os.makedirs(out_path)
    os.system('python GPSnet_module_generation.py %s %s %s &' % (input_path+f,INPUT_PPI_FILE,out_path))

