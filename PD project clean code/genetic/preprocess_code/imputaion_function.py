import pandas as pd
import numpy as np

###### Impute missing
def impute_column(c,df):#column, df
    #c=column_name
    imputed_value = df[['PATNO', 'EVENT_ID', c]]
    pop_median = np.nanmedian(df[c].astype(float))  # population median

    missing_samples = list(set(df.loc[df[c].isna()].PATNO))  # samples that miss the variable
    for s in missing_samples:  # for each sample
        one_sample_df = df.loc[df.PATNO == s][['EVENT_ID', 'INFODT', c]]
        one_sample_df['date'] = pd.to_datetime(one_sample_df['INFODT'],format='%m/%Y')
        one_sample_df_no_missing = one_sample_df.dropna()

        if one_sample_df_no_missing.shape[0] > 0:  # has at least one value
            missing_visits = one_sample_df.loc[one_sample_df[c].isna()][['EVENT_ID', 'date']]
            dates_has_value = one_sample_df_no_missing['date']
            for i in range(missing_visits.shape[0]):
                v = missing_visits.iloc[i, :]['date']
                event = missing_visits.iloc[i, :]['EVENT_ID']
                time_diff = list(abs(dates_has_value - v))
                nearest_ind = time_diff.index(min(time_diff))  # impute by the nearest

                imputed_value.loc[(imputed_value['PATNO'] == s) & (imputed_value['EVENT_ID'] == event), c] = \
                one_sample_df_no_missing.iloc[nearest_ind, :][c]
        else:  # no value for the variable in the sample
            imputed_value.loc[imputed_value.PATNO == s, c] = pop_median  # impute by population median
    return list(imputed_value[c])
    #df2[c + '_impute'] = imputed_value[c]
