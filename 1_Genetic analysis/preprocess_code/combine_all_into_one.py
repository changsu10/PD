##### combine data for LMM
##### results in 3 files: combined_data_for_lme.csv; combined_GWAS_for_lme.csv; combined_nonGWAS_for_lme.csv
import pandas as pd
import numpy as np
from extraction_function import extract_primary_diag
from cal_additional_scores import *
import scipy.stats as stats
from collections import Counter
####
all_clinical_file='../preprocessed_data/longitudinal_all_samples_remove_dupDates.csv'
all_biospecimen_file='../preprocessed_data/biospecimen_all_samples.csv'
all_static_file='../preprocessed_data/static_all_samples.csv'
patient_visit_file='../preprocessed_data/patient_event_date.csv'
ledd_file='../preprocessed_data/ledd.csv'
gwas_name='PD_combined'
gwas_file='../../../amp_pd_manqi/preprocessed_data/PPMI_'+gwas_name+'_genotype.csv'# or PPMI_PD_genotype.csv

patient_visit_df=pd.read_csv(patient_visit_file)
ledd_df=pd.read_csv(ledd_file)
clinical_df=pd.read_csv(all_clinical_file)
biospecimen_df=pd.read_csv(all_biospecimen_file)

static_df=pd.read_csv(all_static_file)
gwas_df=pd.read_csv(gwas_file)
diag_df=extract_primary_diag(path='../../ppmi_data/')

visit_list = ['BL'] + ["V%02d" % i for i in range(1, 21)]
############# filter on static
### select static columns
static_keep_dict={'PATNO':'PATNO',
             'SXDT':'symptom_date','PDDXDT':'diagnosis_date',
             'ANYFAMPD':'is_family_pd','COHORT':'COHORT',
             'ENROLL_DATE':'ENROLL_DATE','ENROLL_STATUS':'ENROLL_STATUS','ENROLL_AGE':'ENROLL_AGE',
             'BIRTHDT':'BIRTHDT','SEX':'SEX','HANDED':'HANDED','EDUCYRS':'EDUCYRS', 'RAWHITE':'is_white'}
static_df=static_df[list(static_keep_dict.keys())]
static_df.rename(columns=static_keep_dict, inplace=True)
print('original data: '+str(static_df.shape[0])+' samples')
## remove withdraw patients
static_df=static_df[static_df.ENROLL_STATUS.isin(['Enrolled','Complete','Screened','Baseline'])]#2381
static_df.drop(columns=['ENROLL_STATUS'],inplace=True)
print('filter on enroll status: '+str(static_df.shape[0])+' samples')
## use PD samples only
static_df=static_df[static_df.COHORT==1]
print('keep PD cohort only: '+str(static_df.shape[0])+' samples')
## use PD diagnosis only
diag1=diag_df.loc[diag_df.PATNO.isin(static_df.PATNO)]
id1=list(set(list(diag1.loc[diag1.PRIMDIAG==1]['PATNO'])))#ids whose diagnosis is PD
static_df=static_df.loc[static_df.PATNO.isin(id1)]
print('keep PD diagnosis only: '+str(static_df.shape[0])+' samples')

### reformat GWAS - add PATNO
patno=[int(i.split('PP-')[-1]) for i in gwas_df['Unnamed: 0']]
gwas_df['PATNO']=patno
gwas_df.drop(columns=['Unnamed: 0'],inplace=True)
## merge wit genetic data, keep samples that have genetic data
df1 = static_df.merge(gwas_df, on=['PATNO'], how='inner')
print('keep samples if has genetic data: '+str(df1.shape[0])+' samples')

## filtered sample ids
sample_ids=list(df1.PATNO)

############ Merge longitudinal
patient_visit_df=patient_visit_df.loc[patient_visit_df.EVENT_ID.isin(visit_list)]
patient_visit_df=patient_visit_df.loc[patient_visit_df.PATNO.isin(sample_ids)]

ledd_df=ledd_df.loc[ledd_df.EVENT_ID.isin(visit_list)]
ledd_df=ledd_df.loc[ledd_df.PATNO.isin(sample_ids)]

biospecimen_df=biospecimen_df.loc[biospecimen_df.EVENT_ID.isin(visit_list)]
biospecimen_df=biospecimen_df.loc[biospecimen_df.PATNO.isin(sample_ids)]

clinical_df=clinical_df.loc[clinical_df.PATNO.isin(sample_ids)]
clinical_df=clinical_df.loc[clinical_df.EVENT_ID.isin(visit_list)]

### Impute moca BL value
# Almost all moca are NA at BL, which makes further normalization infeasible. Here, we imput moca BL by moca SC 
all_clinical_file0='../preprocessed_data/longitudinal_all_samples.csv'#data contains SC
all_clinical_df=pd.read_csv(all_clinical_file0)
all_clinical_df=all_clinical_df.loc[all_clinical_df.PATNO.isin(sample_ids)]
moca_sc_df=all_clinical_df.loc[all_clinical_df.EVENT_ID=='SC'][['PATNO','EVENT_ID','INFODT','MCATOT']]#moca SC df

moca_noSC_df=clinical_df[['PATNO','EVENT_ID','INFODT','MCATOT']]#non SC moca df
moca_sub_df=pd.concat([moca_noSC_df,moca_sc_df])# all moca df

# impute moca BL by SC
for index, row in moca_noSC_df.iterrows():
    patno = row['PATNO']
    event = row['EVENT_ID']
    if event == 'BL':
        if pd.isnull(row['MCATOT']):  # is BL is na
            ## impute by SC
            # bl_date = row['INFODT']
            one_sample = moca_sub_df.loc[moca_sub_df.PATNO == patno]
            if 'SC' in list(one_sample.EVENT_ID):#if has SC
                imputed_values = one_sample.loc[one_sample.EVENT_ID=='SC']['MCATOT']
                imputed_values.dropna(inplace=True)#remove na
                imputed_values = list(set(imputed_values))#remove duplicate
                if len(imputed_values)==1:
                    imputed_value=imputed_values[0]
                elif len(imputed_values)>1:
                    print('multiple value in SC')
                    imputed_value = np.nan
                else:
                    imputed_value = np.nan
            else:
                imputed_value=np.nan
            # one_sample_hasValue = one_sample.dropna()
            #
            # all_dates = list(one_sample_hasValue['INFODT'])
            # all_dates = pd.to_datetime(all_dates, format='%m/%Y')
            #
            # time_diff = list(abs(all_dates - pd.to_datetime(bl_date, format='%m/%Y')))
            # nearest_ind = time_diff.index(min(time_diff))  # nearest
            # imputed_value = one_sample_hasValue.iloc[nearest_ind, :]['MCATOT']
            moca_noSC_df.at[index, 'MCATOT'] = imputed_value
        else:
            pass  # no need to impute
    else:
        pass

clinical_df['MCATOT_imputeBL']=moca_noSC_df['MCATOT']

### calculate additional scores
## calculate tremor/pigd
tremor_items_columns=['NP2TRMR','NP3PTRMR','NP3PTRML','NP3KTRMR','NP3KTRML','NP3RTARU','NP3RTALU','NP3RTARL','NP3RTALL','NP3RTALJ','NP3RTCON']
pigd_items_columns=['NP2WALK','NP2FREZ','NP3GAIT','NP3FRZGT','NP3PSTBL']
tremor_scores,pigd_scores,ratio,phenotype=cal_tremor(clinical_df,tremor_items_columns,pigd_items_columns)
# append to df
clinical_df['tremor_scores']=tremor_scores
clinical_df['pigd_scores']=pigd_scores
#clinical_df['tremor_pigd_category']=phenotype

## calculate GCO
gco_items=['NP1RTOT','NP1PTOT','NP2PTOT','NP3TOT','MSEADLG','MCATOT']
gco=cal_gco(clinical_df,gco_items)
clinical_df['gco']=gco

### select clinical total scores
clinical_df['updrs1']=clinical_df['NP1RTOT']+clinical_df['NP1PTOT']
clinical_keep_dict={'NP2PTOT':'updrs2','NP3TOT':'updrs3','NP4TOT':'updrs4','MSEADLG':'schwab',
                    'tremor_scores':'tremor_scores','pigd_scores':'pigd_scores',#'tremor_pigd_category':'tremor_pigd_category',
                    'MCATOT_imputeBL':'moca','JLO_TOTCALC':'benton','DVT_DELAYED_RECALL':'hvlt_delayed_recall','DVT_RECOG_DISC_INDEX':'hvlt_recog_disc_index',
                    'DVT_RETENTION':'hvlt_retention','DVT_TOTAL_RECALL':'hvlt_total_recall','DVS_LNS':'lns','DVT_SFTANIM':'semantic_fluency',
                    'DVT_SDM':'symbol_digit','geriatric_total':'gds','stai_total':'stai','scopa_total':'scopa','rem_total':'rem',
                    'ess_total':'ess','quip_total':'quip','updrs1':'updrs1','gco':'gco'}#'TOTAL_CORRECT':'upsit',
clinical_df=clinical_df[['PATNO','EVENT_ID','INFODT']+list(clinical_keep_dict.keys())]
clinical_df.rename(columns=clinical_keep_dict, inplace=True)

### select biomarker
biospecimen_df=biospecimen_df[['PATNO','EVENT_ID']+['alpha_syn', 'total_tau', 'abeta_42', 'p_tau181p']]

### merge longitudinal
pd_new = patient_visit_df.merge(clinical_df, on=['PATNO','EVENT_ID','INFODT'], how='outer')
pd_new = pd_new.merge(ledd_df, on=['PATNO','EVENT_ID','INFODT'], how='outer')
pd_new = pd_new.merge(biospecimen_df, on=['PATNO','EVENT_ID'], how='outer')

## remove INFODT is null
df2 = pd_new.loc[pd_new['INFODT'].notna()]
print('remove visits without date:'+str(len(set(df2.PATNO))))
df2.set_index(['PATNO','EVENT_ID','INFODT'], inplace=True)

###### normalize
### 1. convert raw scores to the percentage of the maximum for that scale
### 2. normalized by population baseline mean and std

## step 1: scale
# trait maximum value
col_max=df2.max(axis=0)
# reverse moca and semantic_fluency
df2['moca']=30-df2['moca']
df2['semantic_fluency']=col_max['semantic_fluency']-df2['semantic_fluency']
df2_impute=df2/col_max
## step2: normalize
df2_impute_bl=df2_impute.loc[(slice(None), 'BL'), :]
bl_mean=df2_impute_bl.mean(axis=0)
bl_std=df2_impute_bl.std(axis=0)
df2_norm=(df2_impute-bl_mean)/bl_std
df2_norm=df2_norm.reset_index()

###### Merge static and longitudinal
df=df2_norm.merge(df1,on=['PATNO'],how='inner')

###### add additional covariated and normalize by z-score
## PD_duration, visit-diagnosis_date, month
pd_duration_days=pd.to_datetime(df.INFODT,format='%m/%Y')-pd.to_datetime(df.diagnosis_date,format='%m/%Y')
pd_duration = stats.zscore(pd_duration_days/np.timedelta64(1, 'M'))
pd_duration_norm = list(pd_duration)
df['pd_duration_norm']=pd_duration_norm

## PD symptom duration, visit-symptom_date, month
pd_duration_days=pd.to_datetime(df.INFODT,format='%m/%Y')-pd.to_datetime(df.symptom_date,format='%m/%Y')
symtom_duration=list(stats.zscore(pd_duration_days/np.timedelta64(1, 'M')))
df['symtom_duration_norm']=symtom_duration

## time, visit-baseline, year
time=[]
for i in range(df.shape[0]):
    patno=df.iloc[i,:]['PATNO']
    event=df.iloc[i,:]['EVENT_ID']
    date=df.iloc[i,:]['INFODT']
    bl=df.loc[(df.PATNO==patno) & (df.EVENT_ID=='BL')]['INFODT']
    diff=pd.to_datetime(date,format='%m/%Y')-pd.to_datetime(bl,format='%m/%Y')#event date - baseline date
    diff=diff/np.timedelta64(1, 'Y')
    time.append(diff.values[0])
df['time']=time

#age at onset, diagnosis date - birthdate, year, and z-score normalize
age_days=pd.to_datetime(df.diagnosis_date,format='%m/%Y')-pd.to_datetime(df.BIRTHDT,format='%m/%Y')
age=list(stats.zscore(age_days/np.timedelta64(1, 'Y')))
df['age_onset_norm']=age

# norm LEDD
ledd_norm=list(stats.zscore(df.LEDD))
df['ledd_norm']=ledd_norm

## append baseline values as a new column
value_columns=['updrs1','updrs2', 'updrs3', 'updrs4', 'schwab','tremor_scores', 'pigd_scores', 
               'moca', 'benton', 'hvlt_delayed_recall',
               'hvlt_recog_disc_index', 'hvlt_retention', 'hvlt_total_recall', 'lns',
               'semantic_fluency', 'symbol_digit', 
               'gds', 'stai', 
               'scopa', 
               'rem','ess', 
               'quip', 'gco', 'LEDD', 'alpha_syn', 'total_tau','abeta_42', 'p_tau181p']#trait columns
for value_col in value_columns:
    bl_df = df[df['EVENT_ID'] == 'BL'][['PATNO', value_col]].rename(columns={value_col: f'{value_col}.BL'})
    df = df.merge(bl_df, on='PATNO', how='left')


#### save all
df.to_csv('../preprocessed_data/'+gwas_name+'_combined_data_for_lme.csv',index=False)
## save separately
gwas_columns=[i for i in df.columns if (i.startswith('rs')) or (i.startswith('chr'))]
df_gwas=df[gwas_columns]
df_gwas.to_csv('../preprocessed_data/'+gwas_name+'_combined_GWAS_for_lme.csv',index=False)
other_columns=[i for i in df.columns if not (i.startswith('rs') or (i.startswith('chr')))]
df_others=df[other_columns]
df_others.to_csv('../preprocessed_data/'+gwas_name+'_combined_nonGWAS_for_lme.csv',index=False)

print('All done!')