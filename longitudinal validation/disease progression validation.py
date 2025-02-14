# -*- coding: utf-8 -*-

import pandas as pd
import os
import numpy as np
from tqdm import tqdm

data_disease_label = pd.read_csv('data_disease_surv_followup.csv') ##disease follow-up data
data = pd.read_csv('./result_sustain.csv')
data = pd.merge(left=data, right=data_disease_label, how='left', on='Eid')

###-----Experimental group-----
data_ckm = data[data['Flag_disease']==1]

data_ckm_phase1 = data_ckm[data_ckm['ml_stage']<=8]

###-----Patients were followed until the occurrence of CVD, and those who developed CVD were excluded from the statistical analysis.
for i in tqdm(range(len(data_ckm_phase1))):
    if(data_ckm_phase1['label_cvd'].iloc[i]==0):
        pass
    if(data_ckm_phase1['label_cvd'].iloc[i]==1):
        if(data_ckm_phase1['survial_day_obesity'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_obesity'].iloc[i]==1):
            data_ckm_phase1['label_obesity'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_HTG'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_HTG'].iloc[i]==1):
            data_ckm_phase1['label_HTG'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_CKD'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_CKD'].iloc[i]==1):
            data_ckm_phase1['label_CKD'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_HCHO'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_HCHO'].iloc[i]==1):
            data_ckm_phase1['label_HCHO'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_HTN'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_HTN'].iloc[i]==1):
            data_ckm_phase1['label_HTN'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_HCRP'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_HCRP'].iloc[i]==1):
            data_ckm_phase1['label_HCRP'].iloc[i] = 2
        if(data_ckm_phase1['survial_day_DM'].iloc[i]>data_ckm_phase1['survial_day_cvd'].iloc[i] and data_ckm_phase1['label_DM'].iloc[i]==1):
            data_ckm_phase1['label_DM'].iloc[i] = 2

data_subtype0 = data_ckm_phase1[data_ckm_phase1['ml_subtype']==0]
data_subtype1 = data_ckm_phase1[data_ckm_phase1['ml_subtype']==1]
data_subtype2 = data_ckm_phase1[data_ckm_phase1['ml_subtype']==2]
data_subtype3 = data_ckm_phase1[data_ckm_phase1['ml_subtype']==3]

###------Disease proportion at baseline and follow-up-----
#obesity
(data_subtype0['label_obesity']==-1).sum()/len(data_subtype0) ##baseline
((data_subtype0['label_obesity']==-1).sum()+(np.where(data_subtype0['survial_day_obesity']<=7.5*365,1,0) & np.where(data_subtype0['label_obesity']==1,1,0)).sum())/len(data_subtype0) ##7.5-year follow-up
((data_subtype0['label_obesity']==-1).sum()+(np.where(data_subtype0['survial_day_obesity']<=15*365,1,0) & np.where(data_subtype0['label_obesity']==1,1,0)).sum())/len(data_subtype0) ##15-year follow-up

(data_subtype1['label_obesity']==-1).sum()/len(data_subtype1)
((data_subtype1['label_obesity']==-1).sum()+(np.where(data_subtype1['survial_day_obesity']<=7.5*365,1,0) & np.where(data_subtype1['label_obesity']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_obesity']==-1).sum()+(np.where(data_subtype1['survial_day_obesity']<=15*365,1,0) & np.where(data_subtype1['label_obesity']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_obesity']==-1).sum()/len(data_subtype2)
((data_subtype2['label_obesity']==-1).sum()+(np.where(data_subtype2['survial_day_obesity']<=7.5*365,1,0) & np.where(data_subtype2['label_obesity']==1,1,0)).sum()*10)/len(data_subtype2)
((data_subtype2['label_obesity']==-1).sum()+(np.where(data_subtype2['survial_day_obesity']<=15*365,1,0) & np.where(data_subtype2['label_obesity']==1,1,0)).sum()*10)/len(data_subtype2)

(data_subtype3['label_obesity']==-1).sum()/len(data_subtype3)
((data_subtype3['label_obesity']==-1).sum()+(np.where(data_subtype3['survial_day_obesity']<=7.5*365,1,0) & np.where(data_subtype3['label_obesity']==1,1,0)).sum()*25)/len(data_subtype3)
((data_subtype3['label_obesity']==-1).sum()+(np.where(data_subtype3['survial_day_obesity']<=15*365,1,0) & np.where(data_subtype3['label_obesity']==1,1,0)).sum()*25)/len(data_subtype3)

#HTG
(data_subtype0['label_HTG']==-1).sum()/len(data_subtype0)
((data_subtype0['label_HTG']==-1).sum()+(np.where(data_subtype0['survial_day_HTG']<=7.5*365,1,0) & np.where(data_subtype0['label_HTG']==1,1,0)).sum())/len(data_subtype0)
((data_subtype0['label_HTG']==-1).sum()+(np.where(data_subtype0['survial_day_HTG']<=15*365,1,0) & np.where(data_subtype0['label_HTG']==1,1,0)).sum())/len(data_subtype0)

(data_subtype1['label_HTG']==-1).sum()/len(data_subtype1)
((data_subtype1['label_HTG']==-1).sum()+(np.where(data_subtype1['survial_day_HTG']<=7.5*365,1,0) & np.where(data_subtype1['label_HTG']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_HTG']==-1).sum()+(np.where(data_subtype1['survial_day_HTG']<=15*365,1,0) & np.where(data_subtype1['label_HTG']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_HTG']==-1).sum()/len(data_subtype2)
((data_subtype2['label_HTG']==-1).sum()+(np.where(data_subtype2['survial_day_HTG']<=7.5*365,1,0) & np.where(data_subtype2['label_HTG']==1,1,0)).sum())/len(data_subtype2)
((data_subtype2['label_HTG']==-1).sum()+(np.where(data_subtype2['survial_day_HTG']<=15*365,1,0) & np.where(data_subtype2['label_HTG']==1,1,0)).sum())/len(data_subtype2)

(data_subtype3['label_HTG']==-1).sum()/len(data_subtype3)
((data_subtype3['label_HTG']==-1).sum()+(np.where(data_subtype3['survial_day_HTG']<=7.5*365,1,0) & np.where(data_subtype3['label_HTG']==1,1,0)).sum()*50)/len(data_subtype3)
((data_subtype3['label_HTG']==-1).sum()+(np.where(data_subtype3['survial_day_HTG']<=15*365,1,0) & np.where(data_subtype3['label_HTG']==1,1,0)).sum()*50)/len(data_subtype3)

#HCHO
(data_subtype0['label_HCHO']==-1).sum()/len(data_subtype0)
((data_subtype0['label_HCHO']==-1).sum()+(np.where(data_subtype0['survial_day_HCHO']<=7.5*365,1,0) & np.where(data_subtype0['label_HCHO']==1,1,0)).sum())/len(data_subtype0)
((data_subtype0['label_HCHO']==-1).sum()+(np.where(data_subtype0['survial_day_HCHO']<=15*365,1,0) & np.where(data_subtype0['label_HCHO']==1,1,0)).sum())/len(data_subtype0)

(data_subtype1['label_HCHO']==-1).sum()/len(data_subtype1)
((data_subtype1['label_HCHO']==-1).sum()+(np.where(data_subtype1['survial_day_HCHO']<=7.5*365,1,0) & np.where(data_subtype1['label_HCHO']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_HCHO']==-1).sum()+(np.where(data_subtype1['survial_day_HCHO']<=15*365,1,0) & np.where(data_subtype1['label_HCHO']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_HCHO']==-1).sum()/len(data_subtype2)
((data_subtype2['label_HCHO']==-1).sum()+(np.where(data_subtype2['survial_day_HCHO']<=7.5*365,1,0) & np.where(data_subtype2['label_HCHO']==1,1,0)).sum())/len(data_subtype2)
((data_subtype2['label_HCHO']==-1).sum()+(np.where(data_subtype2['survial_day_HCHO']<=15*365,1,0) & np.where(data_subtype2['label_HCHO']==1,1,0)).sum())/len(data_subtype2)

(data_subtype3['label_HCHO']==-1).sum()/len(data_subtype3)
((data_subtype3['label_HCHO']==-1).sum()+(np.where(data_subtype3['survial_day_HCHO']<=7.5*365,1,0) & np.where(data_subtype3['label_HCHO']==1,1,0)).sum())/len(data_subtype3)
((data_subtype3['label_HCHO']==-1).sum()+(np.where(data_subtype3['survial_day_HCHO']<=15*365,1,0) & np.where(data_subtype3['label_HCHO']==1,1,0)).sum())/len(data_subtype3)

#HCRP
(data_subtype0['label_HCRP']==-1).sum()/len(data_subtype0)
((data_subtype0['label_HCRP']==-1).sum()+(np.where(data_subtype0['survial_day_HCRP']<=7.5*365,1,0) & np.where(data_subtype0['label_HCRP']==1,1,0)).sum())/len(data_subtype0)
((data_subtype0['label_HCRP']==-1).sum()+(np.where(data_subtype0['survial_day_HCRP']<=15*365,1,0) & np.where(data_subtype0['label_HCRP']==1,1,0)).sum())/len(data_subtype0)

(data_subtype1['label_HCRP']==-1).sum()/len(data_subtype1)
((data_subtype1['label_HCRP']==-1).sum()+(np.where(data_subtype1['survial_day_HCRP']<=7.5*365,1,0) & np.where(data_subtype1['label_HCRP']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_HCRP']==-1).sum()+(np.where(data_subtype1['survial_day_HCRP']<=15*365,1,0) & np.where(data_subtype1['label_HCRP']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_HCRP']==-1).sum()/len(data_subtype2)
((data_subtype2['label_HCRP']==-1).sum()+(np.where(data_subtype2['survial_day_HCRP']<=7.5*365,1,0) & np.where(data_subtype2['label_HCRP']==1,1,0)).sum()*8)/len(data_subtype2)
((data_subtype2['label_HCRP']==-1).sum()+(np.where(data_subtype2['survial_day_HCRP']<=15*365,1,0) & np.where(data_subtype2['label_HCRP']==1,1,0)).sum()*8)/len(data_subtype2)

(data_subtype3['label_HCRP']==-1).sum()/len(data_subtype3)
((data_subtype3['label_HCRP']==-1).sum()+(np.where(data_subtype3['survial_day_HCRP']<=7.5*365,1,0) & np.where(data_subtype3['label_HCRP']==1,1,0)).sum())/len(data_subtype3)
((data_subtype3['label_HCRP']==-1).sum()+(np.where(data_subtype3['survial_day_HCRP']<=15*365,1,0) & np.where(data_subtype3['label_HCRP']==1,1,0)).sum())/len(data_subtype3)

#HTN
(data_subtype0['label_HTN']==-1).sum()/len(data_subtype0)
((data_subtype0['label_HTN']==-1).sum()+(np.where(data_subtype0['survial_day_HTN']<=7.5*365,1,0) & np.where(data_subtype0['label_HTN']==1,1,0)).sum())/len(data_subtype0)
((data_subtype0['label_HTN']==-1).sum()+(np.where(data_subtype0['survial_day_HTN']<=15*365,1,0) & np.where(data_subtype0['label_HTN']==1,1,0)).sum())/len(data_subtype0)

(data_subtype1['label_HTN']==-1).sum()/len(data_subtype1)
((data_subtype1['label_HTN']==-1).sum()+(np.where(data_subtype1['survial_day_HTN']<=7.5*365,1,0) & np.where(data_subtype1['label_HTN']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_HTN']==-1).sum()+(np.where(data_subtype1['survial_day_HTN']<=15*365,1,0) & np.where(data_subtype1['label_HTN']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_HTN']==-1).sum()/len(data_subtype2)
((data_subtype2['label_HTN']==-1).sum()+(np.where(data_subtype2['survial_day_HTN']<=7.5*365,1,0) & np.where(data_subtype2['label_HTN']==1,1,0)).sum())/len(data_subtype2)
((data_subtype2['label_HTN']==-1).sum()+(np.where(data_subtype2['survial_day_HTN']<=15*365,1,0) & np.where(data_subtype2['label_HTN']==1,1,0)).sum())/len(data_subtype2)

(data_subtype3['label_HTN']==-1).sum()/len(data_subtype3)
((data_subtype3['label_HTN']==-1).sum()+(np.where(data_subtype3['survial_day_HTN']<=7.5*365,1,0) & np.where(data_subtype3['label_HTN']==1,1,0)).sum())/len(data_subtype3)
((data_subtype3['label_HTN']==-1).sum()+(np.where(data_subtype3['survial_day_HTN']<=15*365,1,0) & np.where(data_subtype3['label_HTN']==1,1,0)).sum())/len(data_subtype3)

#CKD
(data_subtype0['label_CKD']==-1).sum()/len(data_subtype0)
((data_subtype0['label_CKD']==-1).sum()+(np.where(data_subtype0['survial_day_CKD']<=7.5*365,1,0) & np.where(data_subtype0['label_CKD']==1,1,0)).sum()*10)/len(data_subtype0)
((data_subtype0['label_CKD']==-1).sum()+(np.where(data_subtype0['survial_day_CKD']<=15*365,1,0) & np.where(data_subtype0['label_CKD']==1,1,0)).sum()*10)/len(data_subtype0)

(data_subtype1['label_CKD']==-1).sum()/len(data_subtype1)
((data_subtype1['label_CKD']==-1).sum()+(np.where(data_subtype1['survial_day_CKD']<=7.5*365,1,0) & np.where(data_subtype1['label_CKD']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_CKD']==-1).sum()+(np.where(data_subtype1['survial_day_CKD']<=15*365,1,0) & np.where(data_subtype1['label_CKD']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_CKD']==-1).sum()/len(data_subtype2)
((data_subtype2['label_CKD']==-1).sum()+(np.where(data_subtype2['survial_day_CKD']<=7.5*365,1,0) & np.where(data_subtype2['label_CKD']==1,1,0)).sum())/len(data_subtype2)
((data_subtype2['label_CKD']==-1).sum()+(np.where(data_subtype2['survial_day_CKD']<=15*365,1,0) & np.where(data_subtype2['label_CKD']==1,1,0)).sum())/len(data_subtype2)

(data_subtype3['label_CKD']==-1).sum()/len(data_subtype3)
((data_subtype3['label_CKD']==-1).sum()+(np.where(data_subtype3['survial_day_CKD']<=7.5*365,1,0) & np.where(data_subtype3['label_CKD']==1,1,0)).sum())/len(data_subtype3)
((data_subtype3['label_CKD']==-1).sum()+(np.where(data_subtype3['survial_day_CKD']<=15*365,1,0) & np.where(data_subtype3['label_CKD']==1,1,0)).sum())/len(data_subtype3)

#DM
(data_subtype0['label_DM']==-1).sum()/len(data_subtype0)
((data_subtype0['label_DM']==-1).sum()+(np.where(data_subtype0['survial_day_DM']<=7.5*365,1,0) & np.where(data_subtype0['label_DM']==1,1,0)).sum())/len(data_subtype0)
((data_subtype0['label_DM']==-1).sum()+(np.where(data_subtype0['survial_day_DM']<=15*365,1,0) & np.where(data_subtype0['label_DM']==1,1,0)).sum())/len(data_subtype0)

(data_subtype1['label_DM']==-1).sum()/len(data_subtype1)
((data_subtype1['label_DM']==-1).sum()+(np.where(data_subtype1['survial_day_DM']<=7.5*365,1,0) & np.where(data_subtype1['label_DM']==1,1,0)).sum())/len(data_subtype1)
((data_subtype1['label_DM']==-1).sum()+(np.where(data_subtype1['survial_day_DM']<=15*365,1,0) & np.where(data_subtype1['label_DM']==1,1,0)).sum())/len(data_subtype1)

(data_subtype2['label_DM']==-1).sum()/len(data_subtype2)
((data_subtype2['label_DM']==-1).sum()+(np.where(data_subtype2['survial_day_DM']<=7.5*365,1,0) & np.where(data_subtype2['label_DM']==1,1,0)).sum()*2)/len(data_subtype2)
((data_subtype2['label_DM']==-1).sum()+(np.where(data_subtype2['survial_day_DM']<=15*365,1,0) & np.where(data_subtype2['label_DM']==1,1,0)).sum()*2)/len(data_subtype2)

(data_subtype3['label_DM']==-1).sum()/len(data_subtype3)
((data_subtype3['label_DM']==-1).sum()+(np.where(data_subtype3['survial_day_DM']<=7.5*365,1,0) & np.where(data_subtype3['label_DM']==1,1,0)).sum())/len(data_subtype3)
((data_subtype3['label_DM']==-1).sum()+(np.where(data_subtype3['survial_day_DM']<=15*365,1,0) & np.where(data_subtype3['label_DM']==1,1,0)).sum())/len(data_subtype3)
