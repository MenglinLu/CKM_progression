# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cbook as cbook
import os
import pandas as pd
from simfuncs import *
from kde_ebm.mixture_model import fit_all_gmm_models, fit_all_kde_models
from kde_ebm import plotting
import warnings
warnings.filterwarnings("ignore",category=cbook.mplDeprecation)

from pySuStaIn.MixtureSustain import MixtureSustain, MixtureSustainData
import sklearn.model_selection
import pylab

import statsmodels.formula.api as smf

###-----1. Subtype and stage inference of patients at follow-up-----
###-----Instantiation of the MixtureSustain model-----
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)

os.chdir(dir_path)

np.random.seed(5)
data = pd.read_csv('./sample_data.csv')

data = data.dropna(subset=['Creatinine', 'CysC', 'Phosphate', 'Urate', 'Urea', 'eGFR',
                    'Glucose', 'HbA1c', 'Albumin','BMI','Waist','CRP','LDLc','TG', 'CHOL',
                    'HDLc','SBP','DBP'], how='any')

biomarkers = ['Creatinine', 'CysC', 'Phosphate', 'Urate', 'Urea', 'eGFR', 
                    'Glucose', 'HbA1c', 'Albumin','BMI','Waist','CRP','LDLc','TG',
                    'HDLc','SBP','DBP','CHOL']

N = len(biomarkers)
M = data.shape[0]

BiomarkerNames = biomarkers

use_parallel_startpoints = True

N_startpoints = 25
N_S_max = 8
N_iterations_MCMC = int(1e5) 

SuStaInLabels = BiomarkerNames

sustainType = 'mixture_GMM' 
dataset_name = 'data_v1'
output_folder = os.path.join(dir_path, dataset_name + '_' + sustainType)

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

data_control = data[data['Flag_disease']==0]

mean_control = np.mean(np.array(data_control[biomarkers]),axis=0)
std_control = np.std(np.array(data_control[biomarkers]),axis=0)

data_stan = (np.array(data[biomarkers])-mean_control)/std_control
data_control_stan = (np.array(data_control[biomarkers])-mean_control)/std_control

IS_decreasing = np.mean(data_stan,axis=0)<np.mean(data_control_stan,axis=0)

data[['eGFR','HDLc','Albumin','Phosphate']] = data[['eGFR','HDLc','Albumin','Phosphate']] * -1

zdata = pd.DataFrame(data,copy=True)

for biomarker in biomarkers:
    mod = smf.ols('%s ~ Age + Sex + Smoker1 + Alcohol1'%biomarker,  # confounder adjustment
                  data=data[data.Flag_disease==0] # fit this model *only* to individuals in the control group
                 ).fit() # fit model    
    
    # get the "predicted" values for all subjects based on the control model parameters
    predicted = mod.predict(data[['Age', 'Sex', 'Smoker1', 'Alcohol1', biomarker]]) 
    
    # calculate our zscore: observed - predicted / SD of the control group residuals
    w_score = (data.loc[:,biomarker] - predicted) / mod.resid.std()
    
    # save zscore back into our new (copied) dataframe
    zdata.loc[:,biomarker] = w_score

mixtures = fit_all_gmm_models(np.array(zdata[biomarkers]), np.array(zdata['Flag_disease']))

L_yes = np.zeros(zdata[biomarkers].shape)
L_no = np.zeros(zdata[biomarkers].shape)

for i in range(N):
    if sustainType == "mixture_GMM":
        L_no[:, i], L_yes[:, i] = mixtures[i].pdf(None, np.array(zdata[biomarkers])[:, i])
    elif sustainType   == "mixture_KDE":
        L_no[:, i], L_yes[:, i] = mixtures[i].pdf(zdata[:, i].reshape(-1, 1))

sustain = MixtureSustain(L_yes, L_no, SuStaInLabels, N_startpoints, N_S_max, N_iterations_MCMC, output_folder, dataset_name, use_parallel_startpoints)

###-----Loading the trained optimal Sustain model-----
s = 3
pickle_filename_s = output_folder + '/pickle_files/' + dataset_name + '_subtype' + str(s) + '.pickle'
pk = pd.read_pickle(pickle_filename_s)

sample_sequence_i = pk['samples_sequence']
sample_f_i = pk['samples_f']

###-----Individual subtype and stage inference-----
data_val = pd.read_csv('data_followup.csv') ### the longitudinal validation cohort: classified as Phase 1 at baseline; without CVD during follow-up

data_val[['Phosphate','eGFR','Albumin','HDLc']] = data_val[['Phosphate','eGFR','Albumin','HDLc']] * -1

zdata_val = pd.DataFrame(data_val,copy=True)

for biomarker in biomarkers:
    mod = smf.ols('%s ~ Age + Sex + Smoker1+Alcohol1'%biomarker, 
                  data=data_val[data_val.Flag_disease==0] ##healthy controls
                 ).fit() # fit model    
    
    # get the "predicted" values for all subjects based on the control model parameters
    predicted = mod.predict(data_val[['Age', 'Sex', 'Smoker1', 'Alcohol1', biomarker]]) 
    
    # calculate our zscore: observed - predicted / SD of the control group residuals
    w_score = (data_val.loc[:,biomarker] - predicted) / mod.resid.std()
    
    # save zscore back into our new (copied) dataframe
    zdata_val.loc[:,biomarker] = w_score

L_yes_val = np.zeros(zdata_val[biomarkers].shape)
L_no_val = np.zeros(zdata_val[biomarkers].shape)

for i in range(N):
    L_no_val[:, i], L_yes_val[:, i] = mixtures[i].pdf(None, np.array(zdata_val[biomarkers])[:, i])

numStages_new = L_yes.shape[1]

sustainData_newData = MixtureSustainData(L_yes_val, L_no_val, numStages_new)
ml_subtype,         \
        prob_ml_subtype,    \
        ml_stage,           \
        prob_ml_stage,      \
        prob_subtype,       \
        prob_stage,         \
        prob_subtype_stage          = sustain.subtype_and_stage_individuals(sustainData_newData, sample_sequence_i, sample_f_i, 100)

zdata_val['ml_subtype'] = ml_subtype
zdata_val['prob_ml_subtype'] = prob_ml_subtype
zdata_val['ml_stage'] = ml_stage
zdata_val['prob_ml_stage'] = prob_ml_stage

###-----2. Biomarker statistics-----
data_val = data_val[data_val['Flag_disease']==1]
data_val = pd.merge(left=data_val ,right=zdata_val[['Eid','ml_subtype','ml_stage']])
data_val['phase_ins1'] = -1
data_val.loc[np.where(data_val['ml_stage']>=0,1,0)&
                                     np.where(data_val['ml_stage']<=8,1,0)==1,'phase_ins1'] = 1
data_val.loc[np.where(data_val['ml_stage']>=9,1,0)&
                                     np.where(data_val['ml_stage']<=12,1,0)==1,'phase_ins1'] = 2
data_val.loc[np.where(data_val['ml_stage']>=13,1,0)&
                                     np.where(data_val['ml_stage']<=20,1,0)==1,'phase_ins1'] = 3

biomarkers_baseline = ['Creatinine_ins0', 'CysC_ins0', 'Phosphate_ins0', 'Urate_ins0', 'Urea_ins0', 'eGFR_ins0',
                    'Glucose_ins0', 'HbA1c_ins0', 'Albumin_ins0','BMI_ins0','Waist_ins0','CRP_ins0','LDLc_ins0','TG_ins0', 'CHOL_ins0',
                    'HDLc_ins0','SBP_ins0','DBP_ins0']

###-----Biomarker statistics of validation cohort at baseline-----
data_baseline = data_val[biomarkers_baseline+['ml_subtype']]

###-----Biomarker statistics of patients remained in Phase 1 at follow-up-----
data_ins1_inphase1 = data_val[data_val['phase_ins1']==1]

###-----Biomarker statistics of patients progressed to Phase 2 at follow-up-----
data_ins1_inphase2 = data_val[data_val['phase_ins1']==2]

###-----Biomarker statistics of patients progressed to Phase 3 at follow-up-----
data_ins1_inphase3 = data_val[data_val['phase_ins1']==3]
