### remove duplicated dates for the same event
### extract patient_visit_date
import pandas as pd
import numpy as np

all_clinical_file='../preprocessed_data/longitudinal_all_samples.csv'
output_dir = '../preprocessed_data/'
visit_list = ['BL'] + ["V%02d" % i for i in range(1, 21)]
####
clinical_df=pd.read_csv(all_clinical_file)
clinical_df=clinical_df.loc[clinical_df.EVENT_ID.isin(visit_list)]
## extract events with duplicate dates
a=clinical_df[['PATNO','EVENT_ID','INFODT']]
b=a.drop_duplicates(['PATNO','EVENT_ID'])
c=a.merge(b.drop_duplicates(), on=['PATNO','EVENT_ID','INFODT'], how='left', indicator=True)## patients have duplicate events
rep=c.loc[c['_merge'] == 'left_only']# events with duplicate dates
rep=rep.drop_duplicates(['PATNO','EVENT_ID'])
##
clinical_clean_nodup = []#clean parts with duplicated dates
for i in range(rep.shape[0]):
    patno = rep.iloc[i, :]['PATNO']
    event = rep.iloc[i, :]['EVENT_ID']
    # one patient one event records that has duplicated dates
    r = clinical_df.loc[(clinical_df.PATNO == patno) & (clinical_df.EVENT_ID == event)]

    ## keep non-NA values
    new_value = []
    rr = r.drop(columns=['PATNO', 'EVENT_ID', 'INFODT'])#keep numeric columns
    for c in rr.columns:
        v = list(set(list(rr[c].dropna())))
        if len(v) > 1:  # error if have multiple value for the same variable
            print('!!Error!!Multiple records for the same variable!')
            print('%s %s %s' % (str(patno),event,c))
        elif len(v) == 0:  # if all na
            new_value.append(np.nan)
        else:
            new_value += v

    ## keep the correct date
    ## check if the date is an error, shoule be later then the previous event but eariler than the later event
    dup_dates = list(r.INFODT)
    sub_one_sample = clinical_df.loc[clinical_df.PATNO == patno]
    sub_one_sample = sub_one_sample.reset_index(drop=True)
    # date before and after the event
    ind_dup = sub_one_sample.index[sub_one_sample.EVENT_ID == event].tolist()
    date_before = sub_one_sample.iloc[max(min(ind_dup) - 1, 0), :]['INFODT']  # date of the previous event
    date_after = sub_one_sample.iloc[min(max(ind_dup) + 1, sub_one_sample.shape[0] - 1), :]['INFODT']  # date of the later event

    # remove wrong dates
    for d in dup_dates:
        if (pd.to_datetime(date_before) <= pd.to_datetime(d) <= pd.to_datetime(date_after)):#correct date
            pass
        else:# remove wrong date
            dup_dates.remove(d)
            r = r.drop(r.index[r.INFODT == d].tolist())# remove records of the wrong date

    if len(dup_dates) == 1:  # only one correct dates left
        new_row = [patno, event, dup_dates[0]] + new_value

    elif len(dup_dates) == 0:
        print('no correct date for the event')
    else:  # have more than 1 correct dates
        # keep the date with most non-NA records
        rr = r.drop(columns=['PATNO', 'EVENT_ID', 'INFODT'])
        na_count = rr.isnull().sum(axis=1).tolist()
        new_date = dup_dates[na_count.index(min(na_count))]
        new_row = [patno, event, new_date] + new_value

    clinical_clean_nodup.append(new_row)

clinical_clean_nodup = pd.DataFrame(clinical_clean_nodup)
clinical_clean_nodup.columns = list(clinical_df.columns)

## get sub df with no duplicate dates
b=a.drop_duplicates(['PATNO','EVENT_ID'],keep=False)#remove all duplicates
c=a.merge(b.drop_duplicates(), on=['PATNO','EVENT_ID','INFODT'], how='left', indicator=True)
clinical_df_no_dup=clinical_df.iloc[c.index[c._merge=='both'],:]#sub df with no duplicate dates

clinical_df_new=pd.concat([clinical_df_no_dup,clinical_clean_nodup],ignore_index=True)
clinical_df_new.sort_values(by=['PATNO', 'EVENT_ID', 'INFODT'],inplace=True)
clinical_df_new.to_csv(output_dir+'longitudinal_all_samples_remove_dupDates.csv',index=False)
print('Finish longitudinal files!')

### save visit dates
df_sub=clinical_df_new[['PATNO','EVENT_ID','INFODT']]
df_sub2=df_sub[df_sub.EVENT_ID.isin(visit_list)]
df_sub2.drop_duplicates(['PATNO','EVENT_ID','INFODT'],inplace=True)
df_sub2.to_csv(output_dir+'patient_event_date.csv',index=False)

print('Finish all!')

