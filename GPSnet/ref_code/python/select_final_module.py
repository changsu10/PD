import pandas as pd
import numpy as np
import pickle
import sys
import os

trait=sys.argv[1]
top_gene_count_cutoff = float(sys.argv[2])#0.005
save_final_module_path=sys.argv[3]
save_final_module_with_score_path=sys.argv[4]
summary_df=sys.argv[5]

if not os.path.exists(save_final_module_path):
    os.makedirs(save_final_module_path)

if not os.path.exists(save_final_module_with_score_path):
    os.makedirs(save_final_module_with_score_path)

if not os.path.exists(summary_df):
    with open(summary_df, 'w') as f:
        f.write('trait\tcutoff\tsize\n')
    f.close()

# ---------------------------------------------------------------------------
combined_module = {}
combined_score = {}
for start_index in [0,10000,20000,30000,40000,50000,60000,70000,80000,90000]:
    score=pickle.load(open('../Raw_module/'+trait+'/score_'+str(start_index)+'.pkl','rb'))
    raw_module=pickle.load(open('../Raw_module/'+trait+'/module_'+str(start_index)+'.pkl','rb'))
    combined_module.update(raw_module)
    combined_score.update(score)


top_module_score_cutoff = 5 / 100

# remove modules with less than 10 genes
filtered_raw_module = {k: v for k, v in combined_module.items() if len(v) >= 10}

# selec top 5% modules (module score)
sorted_scores = sorted(combined_score.items(), key=lambda item: item[1], reverse=True)
top_n = int(len(combined_score) * top_module_score_cutoff)

# Get the top 5% keys
top_keys = [k for k, v in sorted_scores[:top_n]]

# Select the corresponding sub-modules
top_raw_module = {k: filtered_raw_module[k] for k in top_keys if k in filtered_raw_module}

module_gene = np.concatenate([top_raw_module[idx] for idx in top_raw_module.keys()])

# Unique genes
umg, counts = np.unique(module_gene, return_counts=True)
umg = pd.DataFrame({'gene': umg, 'count': counts})

# Sort by gene significance (count) in descending order
umg.sort_values(by='count', ascending=False, inplace=True)
umg=umg.reset_index()

# Determine cutoff based on the top modules
cutoff_index = umg[umg['count'] < top_gene_count_cutoff * len(top_raw_module)].index.min()
Final_Module = umg.loc[:cutoff_index, 'gene'] if cutoff_index is not None else umg['gene']
Final_Module_with_score=umg.loc[:cutoff_index, ['gene','count']] if cutoff_index is not None else umg['gene']

# Save
Final_Module.to_csv(save_final_module_path+trait+'_'+str(top_gene_count_cutoff)+'.txt',index=False)
Final_Module_with_score.to_csv(save_final_module_with_score_path+trait+'_'+str(top_gene_count_cutoff)+'.txt',index=False)

f=open(summary_df,'a')
f.write(f'{trait}\t{top_gene_count_cutoff}\t{Final_Module.shape[0]}\n')
f.close()