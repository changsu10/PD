from extraction_function import *
import pandas as pd
import math

path='../../ppmi_data/'
version='_05Sep2023'
output_dir='../preprocessed_data/'

static_datasets_type = {
    # static (one patient should have one record only)
    'pd_start': 'medical',
    'family_history': 'subject_characteristics',
    'status': 'subject_characteristics',
    'screening': 'subject_characteristics',
    'socio': 'subject_characteristics',

}
longitudinal_datasets_type={
    # clinical (longitudinal)
    'updrs1': 'motor',
    'updrs1pq': 'motor',
    'updrs2pq': 'motor',
    'updrs3': 'motor',
    'updrs4': 'motor',
    'schwab': 'motor',

    'aut': 'non-motor',
    'cog_catg': 'non-motor',
    'geriatric': 'non-motor',
    'quip': 'non-motor',
    'stai': 'non-motor',
    'benton': 'non-motor',
    'hopkins_verbal': 'non-motor',
    'letter_seq': 'non-motor',
    'moca': 'non-motor',
    'semantic': 'non-motor',
    'sdm': 'non-motor',
    'upsit': 'non-motor',
    'epworth': 'non-motor',
    'rem': 'non-motor',

    'vital_signs': 'medical',
    'neuro_cranial': 'medical',
    #'primary_diag': 'medical',
}
biospecimen_datasets_type={
    # biomarker (longitudinal, format different than clinical)
    'hemoglobin': 'biospecimen',
    'alpha_syn': 'biospecimen',
    'total_tau': 'biospecimen',
    'abeta_42': 'biospecimen',
    'p_tau181p': 'biospecimen',
    'dna': 'biospecimen',
    'rna': 'biospecimen',
    'plasma': 'biospecimen',
    'serum': 'biospecimen',
    'saa': 'biospecimen',
}

feature_vocabulary = {} # save features
feature_source = {}
feature_list = []


combined_static_df=pd.DataFrame()
combined_longitudinal_df=pd.DataFrame()
combined_biospecimen_df=pd.DataFrame()


def append_feature_dict(variable_list,feature_vocabulary,feature_source,feature_list,dataset,datatype):
    ## append feature vocabulary
    for var in variable_list:
        if var in feature_vocabulary:
            print("!!! Error: feature got same name!")
            print(var, feature_source[var])
            print(dataset)
        else:
            feature_vocabulary[var] = datatype
            feature_source[var] = dataset
            feature_list.append(var)
    return feature_vocabulary, feature_source, feature_list


######### ---static features--- ##################
datasets=static_datasets_type.keys()
# extract dataframes
pd_start = extract_pd_start(path,version)
family_history = extract_family_history(path,version)
status = extract_status(path,version)
screening = extract_screening(path,version)
socio = extract_socio(path,version)

for dataset in datasets:
    datatype = static_datasets_type[dataset]
    dataset_df = eval(dataset).reset_index()
    variable_list = [x for x in list(dataset_df.columns) if x not in ['PATNO']]

    feature_vocabulary, feature_source, feature_list = append_feature_dict(variable_list,feature_vocabulary,feature_source,feature_list,dataset,datatype)

    ### combine dataframe
    if combined_static_df.shape[0] == 0:
        combined_static_df = dataset_df
    else:
        common_cols = [value for value in list(dataset_df.columns) if value in list(combined_static_df.columns)]
        if len(common_cols) > 1:
            print('!! Error: more than one common columns except PATNO')
            print(dataset, common_cols)
        combined_static_df = combined_static_df.merge(dataset_df, on=['PATNO'], how='outer')
        combined_static_df.reset_index().drop_duplicates(inplace=True)

combined_static_df.to_csv(output_dir+'static_all_samples.csv',index=False)

print('Finish static files!')



################# ---longitudinal clinical features---#######################
longitudinal_datasets = longitudinal_datasets_type.keys()
#extract all dataframes
updrs1 = extract_updrs1(path,version)
updrs1pq = extract_updrs1pq(path,version)
updrs2pq = extract_updrs2pq(path,version)
updrs3 = extract_updrs3(path,version)
updrs4 = extract_updrs4(path,version)
schwab = extract_schwab(path,version)
aut = extract_aut(path,version)
cog_catg = extract_cog_catg(path,version)
geriatric = extract_geriatric(path,version)
quip = extract_quip(path,version)
stai = extract_stai(path,version)
benton = extract_benton(path,version)
hopkins_verbal = extract_hopkins_verbal(path,version)
letter_seq = extract_letter_seq(path,version)
moca = extract_moca(path,version)
semantic = extract_semantic(path,version)
sdm = extract_sdm(path,version)
upsit = extract_upsit(path,version)
epworth = extract_epworth(path,version)
rem = extract_rem(path,version)
vital_signs = extract_vital_sign(path,version)
neuro_cranial = extract_neuro_cranial(path,version)

for dataset in longitudinal_datasets:
    datatype = longitudinal_datasets_type[dataset]
    dataset_df = eval(dataset)#.reset_index()
    variable_list = [x for x in list(dataset_df.columns) if x not in ['PATNO', 'EVENT_ID', 'INFODT']]

    feature_vocabulary, feature_source, feature_list = append_feature_dict(variable_list, feature_vocabulary,
                                                                           feature_source, feature_list, dataset,
                                                                           datatype)
    ### combine dataframe
    if combined_longitudinal_df.shape[0] == 0:
        combined_longitudinal_df = dataset_df
    else:
        common_cols = [value for value in list(dataset_df.columns) if value in list(combined_longitudinal_df.columns)]
        if len(common_cols) != 3:
            print('!! Error: more than one common columns except PATNO, EVENT_ID, INFODT')
            print(dataset, common_cols)
        combined_longitudinal_df = combined_longitudinal_df.merge(dataset_df, on=common_cols, how='outer')
        combined_longitudinal_df.reset_index().drop_duplicates(inplace=True)

######## save
combined_longitudinal_df.sort_values(by=['PATNO', 'EVENT_ID', 'INFODT'],inplace=True)
combined_longitudinal_df.to_csv(output_dir+'longitudinal_all_samples.csv',index=False)
print('Finish longitudinal files!')

################ ---longitudinal biospecimen features---#######################
## biospecimen - CSF
_, alpha_syn, total_tau, abeta_42, p_tau181p, _, _, _, _ = extract_biospecimen(path, version)
for dataset in ['alpha_syn', 'total_tau', 'abeta_42', 'p_tau181p']:
    datatype = biospecimen_datasets_type[dataset]
    dataset_df = eval(dataset)
    numberic_testvalue = list(dataset_df['TESTVALUE'])
    spec = [float(str(i).split('<')[-1]) for i in numberic_testvalue]#remove <
    dataset_df['testvalue'] = spec
    dataset_df.rename(columns={'testvalue':dataset}, inplace=True)#rename TESTVALUE to TESTNAME.TESTVALUE
    dataset_df.drop(columns=['TESTNAME','TESTVALUE'],inplace=True)

    variable_list = [x for x in list(dataset_df.columns) if x not in ['PATNO', 'EVENT_ID']]#should only have 1 element
    ## append feature dictionary
    feature_vocabulary, feature_source, feature_list = append_feature_dict(variable_list, feature_vocabulary,
                                                                           feature_source, feature_list, dataset,
                                                                           datatype)

    ## combine dataframe
    if combined_biospecimen_df.shape[0] == 0:
        combined_biospecimen_df = dataset_df
    else:
        common_cols = [value for value in list(dataset_df.columns) if value in list(combined_biospecimen_df.columns)]
        if common_cols != ['PATNO', 'EVENT_ID']:
            print('!! Error: unexppected common columns')
            print(dataset, common_cols)
        combined_biospecimen_df = combined_biospecimen_df.merge(dataset_df, on=common_cols, how='outer')
        combined_biospecimen_df.reset_index().drop_duplicates(inplace=True)

# ### biospecimen - SAA
# dataset = 'saa'
# datatype = biospecimen_datasets_type[dataset]
# dataset_df = extract_saa(path, version)
# variable_list = [x for x in list(dataset_df.columns) if x not in ['PATNO', 'EVENT_ID']]
# rename_col={i:'saa.'+i for i in variable_list}
# dataset_df.rename(columns=rename_col, inplace=True)#rename TESTVALUE to TESTNAME.TESTVALUE
# ## append feature dictionary
# feature_vocabulary, feature_source, feature_list = append_feature_dict(variable_list, feature_vocabulary,
#                                                                        feature_source, feature_list, dataset,
#                                                                         datatype)

## combine dataframe
if combined_biospecimen_df.shape[0] == 0:
    combined_biospecimen_df = dataset_df
else:
    common_cols = [value for value in list(dataset_df.columns) if value in list(combined_biospecimen_df.columns)]
    if common_cols != ['PATNO', 'EVENT_ID']:
        print('!! Error: unexppected common columns')
        print(dataset, common_cols)
    combined_biospecimen_df = combined_biospecimen_df.merge(dataset_df, on=common_cols, how='outer')
    combined_biospecimen_df.reset_index().drop_duplicates(inplace=True)

combined_biospecimen_df.sort_values(by=['PATNO', 'EVENT_ID'],inplace=True)
combined_biospecimen_df.to_csv(output_dir+'biospecimen_all_samples.csv',index=False)
print('Finish biospecimen!')

### write feature dicitonary
with open(output_dir+"feature_dictionary.csv", "w") as wf:
    wf.write("Variable,Var_Type,Source\n")
    for var in feature_list:
        wf.write("%s,%s,%s\n" % (var, feature_vocabulary[var], feature_source[var]))
wf.close()

print('All finished!')





