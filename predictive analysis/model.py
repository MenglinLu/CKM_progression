# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 12:11:36 2024

@author: Lenovo
"""

import pandas as pd
from tqdm import tqdm
from sklearn.model_selection import cross_val_score, KFold
from sklearn.linear_model import LogisticRegression

data = pd.read_csv('data_merge.csv')

data_s = pd.read_csv('./result_sustain.py')[['Eid','ml_subtype','ml_stage']]

data = pd.merge(left=data, right=data_s, how='left',on='Eid')

for task in tqdm(['death']):
    data1 = data[data['label_'+task]>=0]
    data1 = pd.get_dummies(data1, columns=['ml_subtype'])
    featuree_s = ['Sex', 'Age', 'Smoker1', 'CHOL', 'HDLc','SBP',
           'hypertensive_treat_ins0', 'diabete_ever_icd10_ins0','ml_subtype_0.0', 'ml_subtype_1.0', 'ml_subtype_2.0', 'ml_subtype_3.0','ml_stage']

    clf_s1 = LogisticRegression(penalty='l2',max_iter=20000)
    kfold = KFold(n_splits=5, shuffle=True, random_state=42)
    auc_scores = cross_val_score(clf_s1, data1.loc[:,featuree_s], data1.loc[:,'label_'+task], cv=kfold, scoring='roc_auc')

    
