# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 12:11:36 2024

@author: Lenovo
"""

import pandas as pd
import os
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import train_test_split
import ascvd
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from tqdm import tqdm
from sklearn.model_selection import train_test_split
import ascvd
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, average_precision_score, f1_score, roc_auc_score

from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import SMOTE, SMOTENC
from imblearn.under_sampling import RandomUnderSampler

from sklearn.metrics import classification_report
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

def metrics_cls2(y_true, y_pred):
    predictions = [int(i>=0.5) for i in y_pred[:,1]]
    acc = accuracy_score(y_true, predictions)
    precision = precision_score(y_true, predictions)
    recall = recall_score(y_true, predictions)
    f1 = f1_score(y_true, predictions)
    auc = roc_auc_score(y_true, y_pred[:,1])
    auprc = average_precision_score(y_true, y_pred[:,1])
    
    return [acc, precision, recall, f1, auc, auprc]

data = pd.read_csv('E:/CKM1/data/data_merge_mets.csv')

data_s = pd.read_csv('E:/CKM1/GMM_mets_K4/result_s3.csv')[['Eid','ml_subtype','ml_stage']]


featuree = ['Sex', 'Age', 'Smoker1', 'CHOL', 'HDLc','SBP',
       'hypertensive_treat_ins0', 'diabete_ever_icd10_ins0']

data = pd.merge(left=data, right=data_s, how='left',on='Eid')[['Eid','Sex', 'Age', 'Smoker1', 'CHOL', 'HDLc','SBP',
       'hypertensive_treat_ins0', 'diabete_ever_icd10_ins0','ml_subtype','ml_stage']]

data_surv = pd.read_csv('E:/CKM1/data/data_forSurv_all_0407_mortality.csv')

data = pd.merge(left=data,right=data_surv, how='left',on='Eid')

for task in tqdm(['death']):
    data1 = data[data['label_'+task]>=0]
    data1 = data1[(np.where(data1['label_'+task]==1,1,0) | 
                           (np.where(data1['label_'+task]==0,1,0) & np.where(data1['survial_day_'+task]>5*365,1,0))==1)]
    
    data1.loc[(np.where(data1['label_'+task]==1,1,0) & 
              np.where(data1['survial_day_'+task]>12*364,1,0))==1,'label_'+task] = 0
    
    under_sample = RandomUnderSampler(sampling_strategy={0:sum(data1['label_'+task])}, random_state=42)
    data_sample, data_y_sample = under_sample.fit_resample(data1, data1['label_'+task])
           
    data_sample_s = pd.get_dummies(data_sample, columns=['ml_subtype'])
    featuree_s = ['Sex', 'Age', 'Smoker1', 'CHOL', 'HDLc','SBP',
           'hypertensive_treat_ins0', 'diabete_ever_icd10_ins0','ml_subtype_0.0', 'ml_subtype_1.0', 'ml_subtype_2.0', 'ml_subtype_3.0','ml_stage']
       
    train_x, test_x, train_y, test_y = train_test_split(data_sample_s.loc[:,featuree_s+['survial_day_'+task,'label_'+task]], data_sample_s['label_'+task], 
                                                            test_size=0.2, stratify=data_sample_s['label_'+task])
        
    clf_s1 = LogisticRegression(penalty='l2',max_iter=20000)
    clf_s1.fit(train_x[featuree_s],train_x['label_'+task])
    
    proba_s1 = clf_s1.predict_proba(test_x[featuree_s])
        
    metric_res = metrics_cls2(test_y, proba_s1)


    