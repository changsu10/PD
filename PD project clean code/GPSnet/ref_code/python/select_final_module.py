import pandas as pd
import numpy as np
import math
import pickle
import sys
### select top modules
#alpha=0.5
trait=sys.argv[1]#'moca'
result_path=sys.argv[2]#'../../FMNW_interpolation/python/Raw_Modules/pval0.05/'
save_path=sys.argv[3]#'../../FMNW_interpolation/python/Final_Modules/pval0.05/'
cutoff=sys.argv[4]#0.005
cutoff=float(cutoff)

## read data
score_dict=pickle.load(open(result_path+trait+'_score.pkl','rb'))
module_dict=pickle.load(open(result_path+trait+'_module.pkl','rb'))
########### test
#from scipy.io import loadmat
# ## assume we have module.pkl and score.pkl
# mat = loadmat('Raw_Module_moca_FamaMacBeth_NeweyWest_50.mat')
# Module = mat['Module']
# Score = mat['Score']
# module_dict = {idx: subarray[0].flatten().tolist() for idx, subarray in enumerate(Module)}
# score_dict = {i: row.tolist() for i, row in enumerate(Score)}
# with open('data.pickle', 'wb') as f:
#     pickle.dump(mat, f)

Score = pd.DataFrame.from_dict(score_dict, orient='index', columns=['seed', 'score', 'size'])
Score['ID'] = range(0, len(Score))
# Remove modules with less than 10 in the third column of Score
Score = Score[Score['size'] >= 10]
Score.sort_values(by='score', ascending=False, inplace=True)

# Select the top 5% of modules
top_5_percent = Score.head(math.ceil(len(Score) * 0.05))
module_gene = np.concatenate([module_dict[idx] for idx in top_5_percent['ID']])

# Unique genes
umg, counts = np.unique(module_gene, return_counts=True)
umg = pd.DataFrame({'gene': umg, 'count': counts})

# Sort by gene significance (count) in descending order
umg.sort_values(by='count', ascending=False, inplace=True)
umg=umg.reset_index()
# Determine cutoff based on the top modules
cutoff_index = umg[umg['count'] < cutoff * len(top_5_percent)].index.min()
Final_Module = umg.loc[:cutoff_index, 'gene'] if cutoff_index is not None else umg['gene']

Final_Module.to_csv(save_path+trait+'_'+str(cutoff)+'.txt',index=False)

