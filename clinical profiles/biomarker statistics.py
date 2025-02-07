# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 17:09:48 2024

@author: Lenovo
"""


import pandas as pd
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler

data_biomarker = pd.read_csv('data_biomarker_raw.csv')  ###The original biomarker data without z-score transformation
data = pd.read_csv('./result_sustain.csv')
data = pd.merge(left=data[['Eid','ml_subtype','ml_stage','Flag_disease']], right=data_biomarker, how='left', on='Eid')

###-----Experimental group------
result_s3 = data[data['Flag_disease']==1]

biomarker = ['Creatinine','CysC', 'Phosphate', 'Urea', 'eGFR', 'Urate', 
                    'Glucose', 'HbA1c', 'Albumin','CRP','CHOL','LDLc','HDLc','TG',
                    'BMI','Waist', 'SBP','DBP']

biomarker_raw = biomarker

###-----Normalization-----
stan = MinMaxScaler()
result_s3_standard = pd.DataFrame(stan.fit_transform(result_s3[biomarker_raw]),columns=biomarker)
result_s3_standard['ml_subtype'] = list(result_s3['ml_subtype'])
result_s3_standard['ml_stage'] = list(result_s3['ml_stage'])

result_s3_standard['phase'] = -1

result_s3_standard.loc[np.where(result_s3_standard['ml_stage']>=0,1,0)&
                                     np.where(result_s3_standard['ml_stage']<=8,1,0)==1,'phase'] = 1
result_s3_standard.loc[np.where(result_s3_standard['ml_stage']>=9,1,0)&
                                     np.where(result_s3_standard['ml_stage']<=12,1,0)==1,'phase'] = 2
result_s3_standard.loc[np.where(result_s3_standard['ml_stage']>=13,1,0)&
                                     np.where(result_s3_standard['ml_stage']<=20,1,0)==1,'phase'] = 3

group_mean = result_s3_standard.groupby(['ml_subtype','phase']).mean()

group_mean.to_csv('result_s3_biomarker_mean.csv')

