# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 17:09:48 2024

@author: Lenovo
"""


import pandas as pd
import os
import numpy as np
from tqdm import tqdm

data_disease_label = pd.read_csv('data_disease_surv.csv')
data = pd.read_csv('./result_sustain.csv')
data = pd.merge(left=data, right=data_disease_label, how='left', on='Eid')

###Experimental group
data_ckm = data[data['Flag_disease']==1]

disease_li = ['label_obesity', 'label_HCHO','label_HTN','label_DM','label_HTG','label_HCRP','label_CKD']

data_subtype0 = data_ckm[data_ckm['ml_subtype']==0]
data_subtype1 = data_ckm[data_ckm['ml_subtype']==1]
data_subtype2 = data_ckm[data_ckm['ml_subtype']==2]
data_subtype3 = data_ckm[data_ckm['ml_subtype']==3]

data_subtype0_phase1 = data_subtype0.loc[np.where(data_subtype0['ml_stage']>=0,1,0)&
                                     np.where(data_subtype0['ml_stage']<=8,1,0)==1,:]
data_subtype0_phase2 = data_subtype0.loc[np.where(data_subtype0['ml_stage']>=9,1,0)&
                                     np.where(data_subtype0['ml_stage']<=12,1,0)==1,:]
data_subtype0_phase3 = data_subtype0.loc[np.where(data_subtype0['ml_stage']>=13,1,0)&
                                     np.where(data_subtype0['ml_stage']<=20,1,0)==1,:]

data_subtype1_phase1 = data_subtype1.loc[np.where(data_subtype1['ml_stage']>=0,1,0)&
                                     np.where(data_subtype1['ml_stage']<=8,1,0)==1,:]
data_subtype1_phase2 = data_subtype1.loc[np.where(data_subtype1['ml_stage']>=9,1,0)&
                                     np.where(data_subtype1['ml_stage']<=12,1,0)==1,:]
data_subtype1_phase3 = data_subtype1.loc[np.where(data_subtype1['ml_stage']>=13,1,0)&
                                     np.where(data_subtype1['ml_stage']<=20,1,0)==1,:]

data_subtype2_phase1 = data_subtype2.loc[np.where(data_subtype2['ml_stage']>=0,1,0)&
                                     np.where(data_subtype2['ml_stage']<=8,1,0)==1,:]
data_subtype2_phase2 = data_subtype2.loc[np.where(data_subtype2['ml_stage']>=9,1,0)&
                                     np.where(data_subtype2['ml_stage']<=12,1,0)==1,:]
data_subtype2_phase3 = data_subtype2.loc[np.where(data_subtype2['ml_stage']>=13,1,0)&
                                     np.where(data_subtype2['ml_stage']<=20,1,0)==1,:]

data_subtype3_phase1 = data_subtype3.loc[np.where(data_subtype3['ml_stage']>=0,1,0)&
                                     np.where(data_subtype3['ml_stage']<=8,1,0)==1,:]
data_subtype3_phase2 = data_subtype3.loc[np.where(data_subtype3['ml_stage']>=9,1,0)&
                                     np.where(data_subtype3['ml_stage']<=12,1,0)==1,:]
data_subtype3_phase3 = data_subtype3.loc[np.where(data_subtype3['ml_stage']>=13,1,0)&
                                     np.where(data_subtype3['ml_stage']<=20,1,0)==1,:]

###Subtype 1
res_dis_i = []
for dis_i in tqdm(disease_li):
    dis_s1_t1 = (data_subtype0_phase1[dis_i]==-1).sum()/len(data_subtype0_phase1)
    dis_s1_t2 = (data_subtype0_phase2[dis_i]==-1).sum()/len(data_subtype0_phase2)
    dis_s1_t3 = (data_subtype0_phase3[dis_i]==-1).sum()/len(data_subtype0_phase3)
    
    
    res_dis_i.append([dis_i,dis_s1_t1,dis_s1_t2,dis_s1_t3])


res_dis_s0_df = pd.DataFrame(res_dis_i,columns=['Disease','t1','t2','t3'])

###Subtype 2
res_dis_i = []
for dis_i in tqdm(disease_li):
    dis_s1_t1 = (data_subtype1_phase1[dis_i]==-1).sum()/len(data_subtype1_phase1)
    dis_s1_t2 = (data_subtype1_phase2[dis_i]==-1).sum()/len(data_subtype1_phase2)
    dis_s1_t3 = (data_subtype1_phase3[dis_i]==-1).sum()/len(data_subtype1_phase3)
    
    
    res_dis_i.append([dis_i,dis_s1_t1,dis_s1_t2,dis_s1_t3])


res_dis_s1_df = pd.DataFrame(res_dis_i,columns=['Disease','t1','t2','t3'])

###Subtype 3
res_dis_i = []
for dis_i in tqdm(disease_li):
    dis_s1_t1 = (data_subtype2_phase1[dis_i]==-1).sum()/len(data_subtype2_phase1)
    dis_s1_t2 = (data_subtype2_phase2[dis_i]==-1).sum()/len(data_subtype2_phase2)
    dis_s1_t3 = (data_subtype2_phase3[dis_i]==-1).sum()/len(data_subtype2_phase3)
    
    
    res_dis_i.append([dis_i,dis_s1_t1,dis_s1_t2,dis_s1_t3])


res_dis_s2_df = pd.DataFrame(res_dis_i,columns=['Disease','t1','t2','t3'])

###Subtype 4
res_dis_i = []
for dis_i in tqdm(disease_li):
    dis_s1_t1 = (data_subtype3_phase1[dis_i]==-1).sum()/len(data_subtype3_phase1)
    dis_s1_t2 = (data_subtype3_phase2[dis_i]==-1).sum()/len(data_subtype3_phase2)
    dis_s1_t3 = (data_subtype3_phase3[dis_i]==-1).sum()/len(data_subtype3_phase3)
    
    
    res_dis_i.append([dis_i,dis_s1_t1,dis_s1_t2,dis_s1_t3])


res_dis_s3_df = pd.DataFrame(res_dis_i,columns=['Disease','t1','t2','t3'])

