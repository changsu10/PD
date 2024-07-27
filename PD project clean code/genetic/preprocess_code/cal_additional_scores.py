### calculate tremor/pigd scores, and GCO scores
import pandas as pd
import numpy as np
from imputaion_function import *
import math
######### calculate Tremor/PIGD phenotype based on updrs2/3
def cal_tremor(df,tremor_items_columns,pigd_items_columns):
    #tremor_items_columns=['NP2TRMR','NP3PTRMR','NP3PTRML','NP3KTRMR','NP3KTRML','NP3RTARU','NP3RTALU','NP3RTARL','NP3RTALL','NP3RTALJ','NP3RTCON']
    #pigd_items_columns=['NP2WALK','NP2FREZ','NP3GAIT','NP3FRZGT','NP3PSTBL']
    #tremor_scores
    tremor_sub=df[['PATNO','EVENT_ID','INFODT']+tremor_items_columns]
    tremor_scores_names=[]
    for c in tremor_items_columns:
        # if tremor_sub[c].isna().sum()>0:# has missing values
        #     new=impute_column(c,tremor_sub)#impute
        #     tremor_sub[c+'_impute']=new
        #     tremor_scores_names.append(c+'_impute')
        # else:
        tremor_scores_names.append(c)
    tremor_scores_sub=tremor_sub[tremor_scores_names]
    tremor_scores=[np.mean(list(tremor_scores_sub.iloc[i,:])) for i in range(tremor_scores_sub.shape[0])]

    #pigd_scores
    pigd_sub=df[['PATNO','EVENT_ID','INFODT']+pigd_items_columns]
    pigd_scores_names = []
    for c in pigd_items_columns:
        # if pigd_sub[c].isna().sum() > 0:
        #     new = impute_column(c,pigd_sub)  # impute
        #     pigd_sub[c + '_impute'] = new
        #     pigd_scores_names.append(c + '_impute')
        # else:
        pigd_scores_names.append(c)
    pigd_scores_sub=pigd_sub[pigd_scores_names]
    pigd_scores=[np.mean(list(pigd_scores_sub.iloc[i,:])) for i in range(pigd_scores_sub.shape[0])]

    # calculate ratio and phenotype
    phenotype=[]
    ratio=[]#tremor/pigd
    for i in range(len(pigd_scores)):
        t=tremor_scores[i]
        p=pigd_scores[i]
        if math.isnan(t) | math.isnan(p):#nan
            ratio.append(np.nan)
            phenotype.append(np.nan)
        else:
            if p==0:
                ratio.append(np.nan)
                if t==0:
                    phenotype.append('IND')
                else:#>0
                    phenotype.append('TD')
            else:
                r=t/p
                ratio.append(r)
                if r >= 1.15:
                    phenotype.append('TD')
                elif r <= 0.9:
                    phenotype.append('PIGD')
                else:
                    phenotype.append('IND')
    return tremor_scores,pigd_scores,ratio,phenotype

############### calculate GCO
#updrs1+updrs2+updrs3+Schwab+moca, 5 vectots,
def cal_gco(df, gco_items):
    #gco_items=['NP1RTOT','NP1PTOT','NP2PTOT','NP3TOT','MSEADLG','MCATOT']
    df_sub=df[['PATNO','EVENT_ID','INFODT']+gco_items]
    gco_score_names=[]
    for c in gco_items:
        if df_sub[c].isna().sum() > 0:
            new = impute_column(c, df_sub)  # impute
            df_sub[c + '_impute'] = new
            gco_score_names.append(c + '_impute')
        else:
            gco_score_names.append(c)

    df_sub=df_sub[['PATNO','EVENT_ID']+gco_score_names]
    df_sub_BL=df_sub.loc[df_sub.EVENT_ID=='BL']

    bl_std=np.std(df_sub_BL[gco_score_names],axis=0)
    bl_mean=np.mean(df_sub_BL[gco_score_names],axis=0)

    gco_items_zscores=[]
    for i in gco_score_names:
        zscore=(df_sub[i]-bl_mean[i])/bl_std[i]
        gco_items_zscores.append(zscore)

    gco_items_zscores=np.array(gco_items_zscores).T
    gco=np.mean(gco_items_zscores,axis=1)#
    return gco
