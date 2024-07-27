import pandas as pd
import numpy as np

def calc_ledd(cm, pn, med_date):
    conmed = cm[cm.PATNO == pn]
    # get the active drugs for patient pn
    active = conmed[(conmed.STARTDT < med_date) & ((conmed.STOPDT >= med_date) | (pd.isnull(conmed.STOPDT)))]
    tot_add_meds = 0
    # check if any of the active meds are a function of levodopa dose
    if (active['LEDD'].str.contains('x') == True).any():
        # get rows that are a function of levodopa dose
        x_spec = active.iloc[np.where(active['LEDD'].str.contains('x') == True)]
        # remove rows from active that are a function of dose
        active = active.drop(x_spec.index.values)
        for i in np.arange(len(x_spec)):
            # pull all active levodopa drugs
            add_meds = active[(active.LEDTRT.str.contains('SINEMET', case=False)) |  # SINEMET
                              (active.LEDTRT.str.contains('CARBIDOPA', case=False)) |  # CARBIDOPA
                              (active.LEDTRT.str.contains('LEVODOPA', case=False)) |  # LEVODOPA
                              (active.LEDTRT.str.contains('RYTARY', case=False)) |  # RYTARY
                              (active.LEDTRT.str.contains('STALEVO', case=False)) |  # STALEVO
                              (active.LEDTRT.str.contains('MADOPAR', case=False)) |  # MADOPAR
                              (active.LEDTRT.str.contains('BENSERAZID', case=False))]  # BENSERAZID
            # pull multiplicative factor
            factor = pd.to_numeric(x_spec.iloc[i].LEDD[-4:])
            # calculate LEDD contribution from levodopa function drugs
            tot_add_meds += pd.to_numeric(add_meds.LEDD).sum() * factor

    return pd.to_numeric(active.LEDD).sum() + tot_add_meds

ledd_df=pd.read_csv('../../ppmi_data/Medical_History/LEDD_Concomitant_Medication_Log_05Sep2023.csv')
patient_visit_df=pd.read_csv('../preprocessed_data/patient_event_date.csv')

ledd_df.STARTDT = pd.to_datetime(ledd_df.STARTDT)
ledd_df.STOPDT = pd.to_datetime(ledd_df.STOPDT)

ledd_processed = pd.DataFrame(columns=['PATNO','EVENT_ID','INFODT','LEDD'])

for idx in range(patient_visit_df.shape[0]):
    pn = patient_visit_df.loc[idx].PATNO
    ei = patient_visit_df.loc[idx].EVENT_ID
    idt = patient_visit_df.loc[idx].INFODT

    ledd = calc_ledd(ledd_df, pn, pd.to_datetime(idt))
    ledd_processed = pd.concat([ledd_processed, pd.DataFrame([pd.Series({'PATNO':pn, 'EVENT_ID':ei, 'INFODT': idt, 'LEDD':ledd})])], ignore_index=True)

ledd_processed=ledd_processed.drop_duplicates()
ledd_processed.to_csv('../preprocessed_data/ledd.csv',index=False)

print('Finish ledd!')