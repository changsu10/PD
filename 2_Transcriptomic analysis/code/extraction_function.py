'''
Functions to extract a dataframe from original file
for longitudinal clinical files:
columns: PATNO, EVENT_ID, INFODT, value1, value2, value2, ...

for longitudinal biospecimen files:
columns: PATNO, EVENT_ID, INFODT, TESTNAME, TESTVALUE

for static files:
columns: value1, value2, ...
index: PATNO
'''

import pandas as pd
import numpy as np
import os
from datetime import datetime

# motor
def extract_updrs1(path='../ppmi_data/',version='_05Sep2023'):
    col_common=["PATNO", "EVENT_ID","INFODT"]#Participant ID, Visit ID, Assessment Date, Date of original data entry
    col_updrs1 = [#"NP1COG", "NP1HALL", "NP1DPRS", "NP1ANXS", "NP1APAT", "NP1DDS",
                  'NP1RTOT']
    updrs1 = pd.read_csv(path+"Motor_Assessments/MDS-UPDRS_Part_I"+version+".csv",
                         index_col=col_common, usecols=col_common+col_updrs1)
    return updrs1.reset_index()

def extract_updrs1pq(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", "INFODT"]
    col_updrs1pq = ["NP1SLPN", "NP1SLPD", "NP1PAIN", "NP1URIN", "NP1CNST", "NP1LTHD", "NP1FATG",
                    'NP1PTOT']
    updrs1pq = pd.read_csv(path+"Motor_Assessments/MDS-UPDRS_Part_I_Patient_Questionnaire"+version+".csv",
                           index_col=col_common, usecols=col_common+col_updrs1pq)
    return updrs1pq.reset_index()

def extract_updrs2pq(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", "INFODT"]
    col_updrs2pq = [#"PATNO", "EVENT_ID",'INFODT',
                    "NP2SPCH", "NP2SALV", "NP2SWAL", "NP2EAT", "NP2DRES", "NP2HYGN","NP2HWRT", "NP2HOBB", "NP2TURN","NP2TRMR", "NP2RISE", "NP2WALK", "NP2FREZ",
                    'NP2PTOT']
    updrs2pq = pd.read_csv(path+"Motor_Assessments/MDS_UPDRS_Part_II__Patient_Questionnaire"+version+".csv",
                           index_col=col_common, usecols=col_common+col_updrs2pq)
    return updrs2pq.reset_index()

def extract_updrs3(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", "INFODT"]
    col_updrs3_temp = [#"PATNO", "EVENT_ID",'INFODT',
                       "PAG_NAME", 'PDMEDYN',#Date exam was administered, PD_MED_USE,
                       "NP3SPCH", "NP3FACXP","NP3RIGN", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL", "NP3FTAPR",
                       "NP3FTAPL", "NP3HMOVR","NP3HMOVL", "NP3PRSPR", "NP3PRSPL", "NP3TTAPR", "NP3TTAPL", "NP3LGAGR",
                       "NP3LGAGL", "NP3RISNG", "NP3GAIT", "NP3FRZGT", "NP3PSTBL", "NP3POSTR", "NP3BRADY", "NP3PTRMR",
                       "NP3PTRML", "NP3KTRMR","NP3KTRML", "NP3RTARU", "NP3RTALU", "NP3RTARL", "NP3RTALL", "NP3RTALJ",
                       "NP3RTCON", "DYSKPRES","DYSKIRAT", "NHY", "HRPOSTMED", "PDSTATE", "PDTRTMNT",
                       'NP3TOT']

    updrs3_temp = pd.read_csv(path+"Motor_Assessments/MDS-UPDRS_Part_III"+version+".csv",
                              index_col=col_common, usecols=col_common+col_updrs3_temp)

    updrs3 = updrs3_temp[updrs3_temp.PAG_NAME == 'NUPDRS3']  # before dose
    #updrs3a = updrs3_temp[updrs3_temp.PAG_NAME == 'NUPDRS3A']  # after dose

    updrs3 = updrs3.reset_index().drop(['PAG_NAME'], axis=1)
    #updrs3a = updrs3a.reset_index().drop(['PAG_NAME'], axis=1)
    return updrs3

def extract_updrs4(path='../ppmi_data/',version='_05Sep2023'):## changed column names
    col_common = ["PATNO", "EVENT_ID", "INFODT"]
    col_updrs4 = ["NP4WDYSK", "NP4DYSKI", "NP4OFF", "NP4FLCTI", "NP4FLCTX", "NP4DYSTN",'NP4WDYSKDEN','NP4WDYSKNUM','NP4WDYSKPCT','NP4DYSTNDEN','NP4OFFNUM','NP4OFFPCT','NP4DYSTNDEN','NP4DYSTNNUM','NP4DYSTNPCT',
                  'NP4TOT']
    updrs4 = pd.read_csv(path+"Motor_Assessments/MDS-UPDRS_Part_IV__Motor_Complications"+version+".csv",
                         index_col=col_common, usecols=col_common+col_updrs4)
    return updrs4.reset_index()

def extract_schwab(path='../ppmi_data/',version='_05Sep2023'):
    col_schwab = ["PATNO", "EVENT_ID", 'INFODT',
                  "MSEADLG"]
    schwab = pd.read_csv(path+"Motor_Assessments/Modified_Schwab___England_Activities_of_Daily_Living"+version+".csv",
                         index_col=["PATNO", "EVENT_ID", 'INFODT'], usecols=col_schwab)
    return schwab.reset_index()

# non-motor
def extract_aut(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", "INFODT"]
    col_aut_1=["SCAU1", "SCAU2", "SCAU3", "SCAU4", "SCAU5", "SCAU6", "SCAU7", "SCAU8", "SCAU9", "SCAU10", "SCAU11",
               "SCAU12", "SCAU13", "SCAU14", "SCAU15", "SCAU16", "SCAU17", "SCAU18", "SCAU19", "SCAU20", "SCAU21"]
    col_aut_2=["SCAU22", "SCAU23","SCAU24", "SCAU25"]

    aut = pd.read_csv(path+"Non-motor_Assessments/SCOPA-AUT"+version+".csv",
                      index_col=col_common, usecols=col_common+col_aut_1+col_aut_2)

    ## In col_aut_1, if value==9, score adds 3. others add value number
    aut1=aut[col_aut_1]
    aut1.replace(9,3,inplace=True)
    aut['total1']=aut1.sum(axis=1)
    ## In col_aut_2, if value==9, score adds 0. Others add value number
    aut2=aut[col_aut_2]
    aut1.replace(9, 0, inplace=True)
    aut['total2'] = aut2.sum(axis=1)
    aut['scopa_total'] = aut['total1']+aut['total2']

    aut = aut[['scopa_total']]

    return aut.reset_index()

def extract_cog_catg(path='../ppmi_data/',version='_05Sep2023'):
    col_cog_catg = ["PATNO", "EVENT_ID", 'INFODT',
                    #"COGDECLN", "FNCDTCOG", "COGSTATE",
                    'COGCAT']
    cog_catg = pd.read_csv(path+"Non-motor_Assessments/Cognitive_Categorization"+version+".csv",
                           index_col=["PATNO", "EVENT_ID"], usecols=col_cog_catg)
    return cog_catg.reset_index()

def extract_geriatric(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_geriatric_pos = ["GDSDROPD", "GDSEMPTY", "GDSBORED", "GDSAFRAD", "GDSHLPLS", "GDSHOME", "GDSMEMRY", "GDSWRTLS", "GDSHOPLS", "GDSBETER"]
    col_geriatric_neg = ["GDSSATIS", "GDSGSPIR", "GDSHAPPY", "GDSALIVE", "GDSENRGY"]

    geriatric = pd.read_csv(path+"Non-motor_Assessments/Geriatric_Depression_Scale__Short_Version_"+version+".csv",
                            index_col=col_common, usecols=col_common+col_geriatric_pos+col_geriatric_neg)

    geriatric["total_pos"] = geriatric[col_geriatric_pos].sum(axis=1)
    sub_neg=1-geriatric[col_geriatric_neg]
    geriatric["total_neg"] = sub_neg.sum(axis=1)

    geriatric["geriatric_total"] = geriatric["total_pos"] + geriatric["total_neg"]

    geriatric = geriatric[["geriatric_total"]]  # drop the rest
    return geriatric.reset_index()

def extract_quip(path='../ppmi_data/',version='_05Sep2023'):
    col_common=["PATNO", "EVENT_ID", 'INFODT']
    col_quip = ["TMGAMBLE", "CNTRLGMB", "TMSEX", "CNTRLSEX", "TMBUY", "CNTRLBUY", "TMEAT", "CNTRLEAT", "TMTORACT", "TMTMTACT", "TMTRWD"]

    quip = pd.read_csv(path+"Non-motor_Assessments/QUIP-Current-Short"+version+".csv",
                       index_col=col_common, usecols=col_common+col_quip)
    quip['quip_total']=quip[col_quip].sum(axis=1)
    quip = quip[['quip_total']]
    return quip.reset_index()

def extract_stai(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_stai = ["STAIAD1", "STAIAD2", "STAIAD3", "STAIAD4", "STAIAD5", "STAIAD6", "STAIAD7", "STAIAD8", "STAIAD9",
                "STAIAD10", "STAIAD11", "STAIAD12", "STAIAD13", "STAIAD14", "STAIAD15", "STAIAD16", "STAIAD17",
                "STAIAD18", "STAIAD19", "STAIAD20", "STAIAD21", "STAIAD22", "STAIAD23", "STAIAD24", "STAIAD25",
                "STAIAD26", "STAIAD27", "STAIAD28", "STAIAD29", "STAIAD30", "STAIAD31", "STAIAD32", "STAIAD33",
                "STAIAD34", "STAIAD35", "STAIAD36", "STAIAD37", "STAIAD38", "STAIAD39", "STAIAD40"]

    col_stai_a_state_pos = ["STAIAD3", "STAIAD4", "STAIAD6", "STAIAD7", "STAIAD9", "STAIAD12", "STAIAD13", "STAIAD14", "STAIAD17", "STAIAD18"]
    col_stai_a_state_neg = ["STAIAD1", "STAIAD2", "STAIAD5", "STAIAD8", "STAIAD10", "STAIAD11", "STAIAD15", "STAIAD16", "STAIAD19", "STAIAD20"]
    col_stai_a_trait_pos = ["STAIAD22", "STAIAD24", "STAIAD25", "STAIAD28", "STAIAD29", "STAIAD31", "STAIAD32","STAIAD35", "STAIAD37", "STAIAD38", "STAIAD40"]
    col_stai_a_trait_neg = ["STAIAD21", "STAIAD23", "STAIAD26", "STAIAD27", "STAIAD30", "STAIAD33", "STAIAD34", "STAIAD36", "STAIAD39"]

    stai = pd.read_csv(path+"Non-motor_Assessments/State-Trait_Anxiety_Inventory"+version+".csv",
                       index_col=col_common, usecols=col_common+col_stai)

    stai["a_state"] = stai[col_stai_a_state_pos].sum(axis=1) + (5 * len(col_stai_a_state_neg) - stai[col_stai_a_state_neg].sum(axis=1))
    stai["a_trait"] = stai[col_stai_a_trait_pos].sum(axis=1) + (5 * len(col_stai_a_trait_neg) - stai[col_stai_a_trait_neg].sum(axis=1))

    stai['stai_total']=stai["a_state"]+stai["a_trait"]
    stai = stai[['stai_total']]

    return stai.reset_index()

def extract_benton(path='../ppmi_data/',version='_05Sep2023'):
    #col_benton = ["PATNO", "EVENT_ID", 'INFODT','JLO_TOTCALC', "JLO_TOTRAW", "DVS_JLO_MSSA", "DVS_JLO_MSSAE"] # "AGE_ASSESS_JLO",
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_benton = ['JLO_TOTCALC']
    benton = pd.read_csv(path+"Non-motor_Assessments/Benton_Judgement_of_Line_Orientation"+version+".csv",
                         index_col=col_common, usecols=col_common+col_benton)
    return benton.reset_index()

def extract_hopkins_verbal(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_hopkins_verbal = ['DVT_DELAYED_RECALL','DVT_RECOG_DISC_INDEX','DVT_RETENTION','DVT_TOTAL_RECALL']
    # col_hopkins_verbal = ["PATNO", "EVENT_ID", 'INFODT', 'DVT_TOTAL_RECALL',
    #                       "HVLTRT1", "HVLTRT2", "HVLTRT3", "HVLTRDLY", "HVLTREC", "HVLTFPRL", "HVLTFPUN"]
    hopkins_verbal = pd.read_csv(path+"Non-motor_Assessments/Hopkins_Verbal_Learning_Test_-_Revised"+version+".csv",
                                index_col=col_common, usecols=col_common+col_hopkins_verbal)
    return hopkins_verbal.reset_index()

def extract_letter_seq(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_letter_seq = [#"PATNO", "EVENT_ID", 'INFODT',
                      'DVS_LNS'#, "LNS_TOTRAW"
                      #"LNS1A", "LNS1B", "LNS1C", "LNS2A","LNS2B", "LNS2C", "LNS3A", "LNS3B", "LNS3C", "LNS4A", "LNS4B",
                      #"LNS4C", "LNS5A", "LNS5B", "LNS5C", "LNS6A", "LNS6B", "LNS6C", "LNS7A", "LNS7B", "LNS7C",
                      ]
    # col_letter_seq_details = ["LNS1A", "LNS1B", "LNS1C", "LNS2A","LNS2B", "LNS2C", "LNS3A", "LNS3B", "LNS3C", "LNS4A",
    #                           "LNS4B", "LNS4C", "LNS5A", "LNS5B", "LNS5C", "LNS6A", "LNS6B", "LNS6C", "LNS7A", "LNS7B", "LNS7C"]

    letter_seq = pd.read_csv(path+"Non-motor_Assessments/Letter_-_Number_Sequencing"+version+".csv",
                             index_col=col_common, usecols=col_common+col_letter_seq)

    #letter_seq["total"] = letter_seq[col_letter_seq_details].sum(axis=1)
    #letter_seq = letter_seq[["INFODT", "total"]]  # letter_seq[["total"]] or letter_seq[["LNS_TOTRAW"]]
    return letter_seq.reset_index()

def extract_moca(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_moca = [#"PATNO", "EVENT_ID", "INFODT",
                #"MCAALTTM", "MCACUBE", "MCACLCKC", "MCACLCKN", "MCACLCKH", "MCALION", "MCARHINO", "MCACAMEL", "MCAFDS","MCABDS", "MCAVIGIL", "MCASER7", "MCASNTNC", "MCAVFNUM", "MCAVF", "MCAABSTR", "MCAREC1", "MCAREC2","MCAREC3", "MCAREC4", "MCAREC5", "MCADATE", "MCAMONTH", "MCAYR", "MCADAY", "MCAPLACE", "MCACITY",
                "MCATOT"]

    # col_moca_visuospatial = [ "MCAALTTM", "MCACUBE", "MCACLCKC", "MCACLCKN", "MCACLCKH"]
    # col_moca_naming = [ "MCALION", "MCARHINO", "MCACAMEL"]
    # col_moca_attention = [ "MCAFDS", "MCABDS", "MCAVIGIL", "MCASER7"]
    # col_moca_language = [ "MCASNTNC", "MCAVF"]
    # col_moca_delayed_recall = [ "MCAREC1", "MCAREC2", "MCAREC3", "MCAREC4", "MCAREC5"]
    # col_moca_orientation = [ "MCADATE", "MCAMONTH", "MCAYR", "MCADAY", "MCAPLACE", "MCACITY"]

    moca = pd.read_csv( path+"Non-motor_Assessments/Montreal_Cognitive_Assessment__MoCA_"+version+".csv",
                        index_col=col_common, usecols=col_common+col_moca)

    # moca["visuospatial"] = moca[col_moca_visuospatial].sum(axis=1)
    # moca["naming"] = moca[col_moca_naming].sum(axis=1)
    # moca["attention"] = moca[col_moca_attention].sum(axis=1)
    # moca["language"] = moca[col_moca_language].sum(axis=1)
    # moca["delayed_recall"] = moca[col_moca_delayed_recall].sum(axis=1)
    #moca = moca[["INFODT", "visuospatial", "naming", "attention", "language", "delayed_recall", "MCAABSTR", "MCAVFNUM", "MCATOT"]]  # drop extra

    return moca.reset_index()

def extract_semantic(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_semantic = ['DVT_SFTANIM']
    semantic = pd.read_csv(path+"Non-motor_Assessments/Modified_Semantic_Fluency"+version+".csv",
                           index_col=col_common, usecols=col_common+col_semantic)
    return semantic.reset_index()

def extract_sdm(path='../ppmi_data/',version='_05Sep2023'):
    col_common = ["PATNO", "EVENT_ID", 'INFODT']
    col_sdm = ["DVT_SDM"]
    sdm = pd.read_csv(path+"Non-motor_Assessments/Symbol_Digit_Modalities_Test"+version+".csv",
                      index_col=col_common, usecols=col_common+col_sdm)
    return sdm.reset_index()

def extract_upsit(path='../ppmi_data/',version='_05Sep2023'):#no Olfactory_UPSIT files anymore, it's in upsit
    col_upsit = ['PATNO', 'EVENT_ID','INFODT', 'SCENT_01_CORRECT',  'SCENT_02_CORRECT', 'SCENT_03_CORRECT',
                  'SCENT_04_CORRECT',  'SCENT_05_CORRECT', 'SCENT_06_CORRECT','SCENT_07_CORRECT',  'SCENT_08_CORRECT',
                  'SCENT_09_CORRECT',  'SCENT_10_CORRECT', 'SCENT_11_CORRECT',  'SCENT_12_CORRECT', 'SCENT_13_CORRECT',
                  'SCENT_14_CORRECT',  'SCENT_15_CORRECT', 'SCENT_16_CORRECT',  'SCENT_17_CORRECT', 'SCENT_18_CORRECT',
                  'SCENT_19_CORRECT', 'SCENT_20_CORRECT', 'SCENT_21_CORRECT', 'SCENT_22_CORRECT', 'SCENT_23_CORRECT',
                  'SCENT_24_CORRECT', 'SCENT_25_CORRECT', 'SCENT_26_CORRECT', 'SCENT_27_CORRECT',  'SCENT_28_CORRECT',
                  'SCENT_29_CORRECT',  'SCENT_30_CORRECT','SCENT_31_CORRECT',  'SCENT_32_CORRECT',  'SCENT_33_CORRECT',
                  'SCENT_34_CORRECT',  'SCENT_35_CORRECT', 'SCENT_36_CORRECT',  'SCENT_37_CORRECT',  'SCENT_38_CORRECT',
                 'SCENT_39_CORRECT', 'SCENT_40_CORRECT', 'TOTAL_CORRECT']
    upsit = pd.read_csv(path+"Non-motor_Assessments/University_of_Pennsylvania_Smell_Identification_Test__UPSIT_"+version+".csv",
                        usecols=col_upsit, index_col=["PATNO", "EVENT_ID",'INFODT'])
    return upsit.reset_index()

def extract_epworth(path='../ppmi_data/',version='_05Sep2023'):
    col_epworth = ["PATNO", "EVENT_ID", 'INFODT',
                   "ESS1", "ESS2", "ESS3", "ESS4", "ESS5", "ESS6", "ESS7", "ESS8"]
    epworth = pd.read_csv(path+"Non-motor_Assessments/Epworth_Sleepiness_Scale"+version+".csv",
                          index_col=["PATNO", "EVENT_ID",'INFODT'], usecols=col_epworth)
    epworth['ess_total']=epworth[["ESS1", "ESS2", "ESS3", "ESS4", "ESS5", "ESS6", "ESS7", "ESS8"]].sum(axis=1)
    epworth=epworth[['ess_total']]
    return epworth.reset_index()

def extract_rem(path='../ppmi_data/',version='_05Sep2023'):
    col_rem = [ "PATNO", "EVENT_ID", 'INFODT',
                "DRMVIVID", "DRMAGRAC", "DRMNOCTB", "SLPLMBMV", "SLPINJUR", "DRMVERBL", "DRMFIGHT", "DRMUMV", "DRMOBJFL",
                "MVAWAKEN", "DRMREMEM", "SLPDSTRB", "STROKE", "HETRA", "PARKISM", "RLS", "NARCLPSY", "DEPRS", "EPILEPSY", "BRNINFM", "CNSOTH" ]

    col_rem_calculate=["DRMVIVID", "DRMAGRAC", "DRMNOCTB", "SLPLMBMV", "SLPINJUR", "DRMVERBL", "DRMFIGHT", "DRMUMV", "DRMOBJFL",
                "MVAWAKEN", "DRMREMEM", "SLPDSTRB", "STROKE", "HETRA", "PARKISM", "RLS", "NARCLPSY", "DEPRS", "EPILEPSY", "BRNINFM", "CNSOTH" ]

    rem = pd.read_csv(path+"Non-motor_Assessments/REM_Sleep_Behavior_Disorder_Questionnaire"+version+".csv",
                      index_col=["PATNO", "EVENT_ID",'INFODT'], usecols=col_rem)

    rem['rem_total']= rem[col_rem_calculate].sum(axis=1)
    rem=rem[['rem_total']]
    return rem.reset_index()

# biospecimen
def extract_biospecimen(path='../ppmi_data/',version='_05Sep2023'):
    col_biospecimen = ["PATNO", "CLINICAL_EVENT", #"RUNDATE",
                       "TYPE", "TESTNAME", "TESTVALUE", "UNITS"]
    biospecimen = pd.read_csv(path+"biospecimen/Current_Biospecimen_Analysis_Results"+version+".csv",
                              index_col=["PATNO"], usecols=col_biospecimen, dtype={'UNITS': str})
    biospecimen.rename(columns={'CLINICAL_EVENT':'EVENT_ID'}, inplace=True)

    csf = biospecimen[(biospecimen["TYPE"] == 'Cerebrospinal Fluid') &
                      ~(biospecimen["TESTVALUE"] == "below detection limit")][["EVENT_ID", "TESTNAME", "TESTVALUE"]]
    hemoglobin = csf[csf["TESTNAME"] == "CSF Hemoglobin"].reset_index().drop_duplicates(["PATNO","EVENT_ID","TESTNAME"])
    alpha_syn = csf[csf["TESTNAME"] == "CSF Alpha-synuclein"].reset_index().drop_duplicates(["PATNO","EVENT_ID","TESTNAME"])
    total_tau = csf[csf["TESTNAME"] == "tTau"].reset_index().drop_duplicates(["PATNO","EVENT_ID","TESTNAME"])
    abeta_42 = csf[csf["TESTNAME"] == "ABeta 1-42"].reset_index().drop_duplicates(["PATNO","EVENT_ID","TESTNAME"])
    p_tau181p = csf[csf["TESTNAME"] == "pTau"].reset_index().drop_duplicates(["PATNO","EVENT_ID","TESTNAME"])
    dna = biospecimen[(biospecimen["TYPE"] == 'DNA')][["EVENT_ID", "TESTNAME", "TESTVALUE"]]
    rna = biospecimen[(biospecimen["TYPE"] == 'RNA')][["EVENT_ID", "TESTNAME", "TESTVALUE"]]
    plasma = biospecimen[(biospecimen["TYPE"] == 'Plasma')][["EVENT_ID", "TESTNAME", "TESTVALUE"]]
    serum = biospecimen[(biospecimen["TYPE"] == 'Serum')][["EVENT_ID", "TESTNAME", "TESTVALUE"]]


    return hemoglobin, alpha_syn, total_tau, abeta_42, p_tau181p, dna, rna, plasma, serum

def extract_saa(path='../ppmi_data/',version='_05Sep2023'):
    col_saa=['PATNO','CLINICAL_EVENT', #'RUNDATE',
             'FmaxRep1', 'FmaxRep2', 'FmaxRep3', 'T50Rep1', 'T50Rep2', 'T50Rep3','TTTRep1', 'TTTRep2', 'TTTRep3',
             'AUCRep1', 'AUCRep2', 'AUCRep3','SLOPERep1', 'SLOPERep2', 'SLOPERep3','InstrumentRep1', 'InstrumentRep2',
             'InstrumentRep3', 'TSmaxRep1','TSmaxRep2', 'TSmaxRep3',
             'SampleVolRep1', 'SampleVolRep2','SampleVolRep3', 'QUALAlgoRep1', 'QUALAlgoRep2', 'QUALAlgoRep3', 'QUALRep1', 'QUALRep2', 'QUALRep3','SynucleinRep1', 'SynucleinRep2', 'SynucleinRep3',#character
             'SLOPEMaxRep1', 'SLOPEMaxRep2', 'SLOPEMaxRep3']
    saa=pd.read_csv(path+"biospecimen/SAA_Biospecimen_Analysis_Results"+version+".csv",index_col=["PATNO"], usecols=col_saa)
    saa.rename(columns={'CLINICAL_EVENT': 'EVENT_ID'}, inplace=True)
    return saa.reset_index()

def extract_vital_sign(path='../ppmi_data/',version='_05Sep2023'):
    col_vital_sign = ["PATNO", "EVENT_ID", 'INFODT',
                      "WGTKG","HTCM"]
    vital_signs = pd.read_csv(path+"Medical_History/Vital_Signs"+version+".csv",
                              index_col=["PATNO", "EVENT_ID"], usecols=col_vital_sign)
    return vital_signs.reset_index()

# Medical-Neurological Exam
def extract_neuro_cranial(path='../ppmi_data/',version='_05Sep2023'):
    col_neuro_cranial = ["PATNO", "EVENT_ID", 'INFODT',
                         'CNRSP']#only found this column
    neuro_cranial = pd.read_csv(path+"Medical_History/Neurological_Exam"+version+".csv",
        index_col=["PATNO", "EVENT_ID"], usecols=col_neuro_cranial)
    return neuro_cranial.reset_index()


# # used to be in enrollment
# # not used in the code?
def extract_primary_diag(path='../ppmi_data/',version='_05Sep2023'):
    col_primary_diag = ["PATNO", 'EVENT_ID', 'INFODT',
                        "PRIMDIAG"]
    primary_diag = pd.read_csv(path+"Medical_History/Primary_Clinical_Diagnosis"+version+".csv",
                               index_col=["PATNO","EVENT_ID"], usecols=col_primary_diag)
    return primary_diag.reset_index()

########## static, one patient only has one record
def extract_pd_start(path='../ppmi_data/',version='_05Sep2023'):#no date
    ## SXMO and SXYEAR are replaced by SXDT
    col_pd_features = ["PATNO",
                       "SXDT", "PDDXDT"] # date of symptom , date of diagnosis
    pd_start = pd.read_csv(path+"Medical_History/PD_Diagnosis_History"+version+".csv", index_col=["PATNO"],
                          usecols=col_pd_features)
    return pd_start


# Subject Characteristics
def extract_family_history(path='../ppmi_data/',version='_05Sep2023'):# no date
    col_family_history = ["PATNO", 'EVENT_ID', 'ORIG_ENTRY',
                          'ANYFAMPD', "BIOMOMPD",  "BIODADPD", "FULSIBPD",  "HAFSIBPD",
                           "MAGPARPD", "PAGPARPD",  "MATAUPD", "PATAUPD",  "KIDSPD"]
    df = pd.read_csv(path+"Subject_Characteristics/Family_History"+version+".csv", usecols=col_family_history )
    ## take BL value, if no BL take the earliest value,
    new_df = []
    all_patients = list(set(list(df.PATNO)))
    all_patients.sort()
    for i in all_patients:
        one_patient = df.loc[df.PATNO == i]
        if one_patient.shape[0] > 1:
            visits = list(one_patient.EVENT_ID)
            if 'BL' in visits: #if has BL
                add = one_patient.loc[one_patient.EVENT_ID == 'BL']
            else:  # choose the latearliestest record
                time = list(one_patient.ORIG_ENTRY)#if use INFODT, there are NAs
                time_new = [datetime.strptime(i, '%M/%Y') for i in time]
                ind = time_new.index(min(time_new))
                add = one_patient.iloc[ind, :]
        else:
            add = one_patient
        add = add.squeeze()
        new_df.append(add)
    df_new = pd.concat(new_df,axis=1,ignore_index = True).T
    df_new['ANYFAMPD'] = df.loc[:, ["BIOMOMPD", "BIODADPD", "FULSIBPD", "HAFSIBPD", "MAGPARPD",
                                    "PAGPARPD", "MATAUPD", "PATAUPD", "KIDSPD"]].sum(axis=1)
    df_new['ANYFAMPD'][df_new['ANYFAMPD']>1]=1
    df_new.drop(columns=['EVENT_ID','ORIG_ENTRY'],inplace=True)
    df_new.set_index('PATNO', inplace=True)
    return df_new

def extract_status(path='../ppmi_data/',version='_05Sep2023'):#no date
    col_status = ["PATNO", "ENROLL_DATE",'COHORT','COHORT_DEFINITION',
                  'ENROLL_AGE','ENROLL_STATUS']#added columns
    status = pd.read_csv(path+"Subject_Characteristics/Participant_Status"+version+".csv",
                         index_col=["PATNO"], usecols=col_status)
    return status

def extract_screening(path='../ppmi_data/',version='_05Sep2023'):#no date
    ## including socio previous handed
    ## APPRDX moved to patient_status COHORT
    ## CURRENT_APPRDX, excluded
    col_screening = ["PATNO",
                     "BIRTHDT", "SEX", "HISPLAT", "HANDED", "RAINDALS", "RAASIAN", "RABLACK", "RAHAWOPI", "RAWHITE", "RANOS"]
    screening = pd.read_csv(path+"Subject_Characteristics/Demographics"+version+".csv",
                index_col=["PATNO"], usecols=col_screening)
    return screening

def extract_socio(path='../ppmi_data/',version='_05Sep2023'):#no date
    col_socio = [ "PATNO", 'EVENT_ID', #'ORIG_ENTRY', 'EVENT_ID',
                  "EDUCYRS"]
    socio = pd.read_csv(path+"Subject_Characteristics/Socio-Economics"+version+".csv", usecols=col_socio)
    ## take BL value, if no BL take SC value,
    all_patient = list(set(socio.PATNO))
    eduyears = []
    for i in all_patient:
        one_record = socio.loc[socio.PATNO == i]
        if one_record.shape[0] > 1:  # more than 2 records
            if 'BL' in list(one_record.EVENT_ID):
                eduyears += list(one_record[one_record.EVENT_ID == 'BL']['EDUCYRS'].values)
            elif 'SC' in list(one_record.EVENT_ID):
                eduyears += list(one_record[one_record.EVENT_ID == 'SC']['EDUCYRS'].values)
            else:#only has 3 possible event in socio: BL, SC, TRANS
                #print(one_record.EVENT_ID)
                print('!!! Error in extract_socio !! Unexpected EVENT_ID')
        else:
            eduyears += list(one_record['EDUCYRS'].values)
    #dic = {'PATNO': all_patient, 'EDUCYRS': eduyears}
    socio_new = pd.DataFrame(eduyears,index=all_patient,columns=['EDUCYRS'])
    socio_new.index.name = 'PATNO'
    socio_new=socio_new.sort_index()
    return socio_new

