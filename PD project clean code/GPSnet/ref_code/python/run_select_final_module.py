import os

cutoff_list=[0.005,0.01,0.05,0.1,0.5]
input_path='../../FMNW_interpolation/python/Raw_Modules/pval0.05/'
save_path='../../FMNW_interpolation/python/Final_Modules/pval0.05/'

file_list=os.listdir(input_path)
file_list=[i for i in file_list if i.endswith('_module.pkl')]
trait_list=[('_').join(i.split('_')[:-1]) for i in file_list]

for t in trait_list:
    for cutoff in cutoff_list:
        cmd='python select_final_module.py %s %s %s %f' % (t,input_path,save_path,cutoff)
        os.system(cmd)
        print('finish '+t)

