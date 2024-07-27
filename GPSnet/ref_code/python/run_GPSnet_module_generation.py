import os
import sys

input_path='../planB/GPSnet_input/per_trait_final/'
DIR_MODULE = '../planB/Raw_Modules/per_trait_final/'
INPUT_PPI_FILE='../ref/ppi.csv'

#if not os.path.exists(INPUT_PPI_FILE):
#    os.makedirs(INPUT_PPI_FILE)

file_list=os.listdir(input_path)
file_list=[i for i in file_list if i.endswith('.csv')]

for f in file_list:
	os.system('python GPSnet_module_generation.py %s %s %s &' % (input_path+f,INPUT_PPI_FILE,DIR_MODULE))

