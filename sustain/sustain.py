import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cbook as cbook
import os
import pandas as pd
from simfuncs import *
from kde_ebm.mixture_model import fit_all_gmm_models, fit_all_kde_models
from kde_ebm import plotting
import warnings
warnings.filterwarnings("ignore",category=cbook.mplDeprecation)
from pySuStaIn.MixtureSustain import MixtureSustain

import statsmodels.formula.api as smf

file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)

os.chdir(dir_path)

np.random.seed(5)

###------Reading data------
data = pd.read_csv('./sample_data.csv')

data = data.dropna(subset=['Creatinine', 'CysC', 'Phosphate', 'Urate', 'Urea', 'eGFR',
                    'Glucose', 'HbA1c', 'Albumin','BMI','Waist','CRP','LDLc','TG', 'CHOL',
                    'HDLc','SBP','DBP'], how='any')

biomarkers = ['Creatinine', 'CysC', 'Phosphate', 'Urate', 'Urea', 'eGFR', 
                    'Glucose', 'HbA1c', 'Albumin','BMI','Waist','CRP','LDLc','TG',
                    'HDLc','SBP','DBP','CHOL']

###------Setting hyperparameters of SuStaIn model------
N = len(biomarkers)
M = data.shape[0]

BiomarkerNames = biomarkers

use_parallel_startpoints = True

N_startpoints = 25
N_S_max = 8
N_iterations_MCMC = int(1e5) 

SuStaInLabels = BiomarkerNames

###------Setting output folder------
sustainType = 'mixture_GMM' 
dataset_name = 'data_v1'
output_folder = os.path.join(dir_path, dataset_name + '_' + sustainType)

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

###------Data preprocessing------
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

###------MixtureSustain model-----
mixtures = fit_all_gmm_models(np.array(zdata[biomarkers]), np.array(zdata['Flag_disease']))

L_yes = np.zeros(zdata[biomarkers].shape)
L_no = np.zeros(zdata[biomarkers].shape)

for i in range(N):
    if sustainType == "mixture_GMM":
        L_no[:, i], L_yes[:, i] = mixtures[i].pdf(None, np.array(zdata[biomarkers])[:, i])
    elif sustainType   == "mixture_KDE":
        L_no[:, i], L_yes[:, i] = mixtures[i].pdf(zdata[:, i].reshape(-1, 1))

sustain = MixtureSustain(L_yes, L_no, SuStaInLabels, N_startpoints, N_S_max, N_iterations_MCMC, output_folder, dataset_name, use_parallel_startpoints)

###-----Running MixtureSustain model------
samples_sequence,   \
    samples_f,          \
    ml_subtype,         \
    prob_ml_subtype,    \
    ml_stage,           \
    prob_ml_stage,      \
    prob_subtype_stage      = sustain.run_sustain_algorithm(plot=True)

###-----Results output based on the optimal number of subtypes-----
###-----After determining the optimal number of subtypes through inference using sustain-cv.py, 
###-----set s = the optimal number of subtypes-1.
s = 3 
pickle_filename_s = output_folder + '/pickle_files/' + dataset_name + '_subtype' + str(s) + '.pickle'
pk = pd.read_pickle(pickle_filename_s)

for variable in ['ml_subtype', # the assigned subtype
                 'prob_ml_subtype', # the probability of the assigned subtype
                 'ml_stage', # the assigned stage 
                 'prob_ml_stage',]: # the probability of the assigned stage
    
    # add SuStaIn output to dataframe
    zdata.loc[:,variable] = pk[variable] 

zdata.to_csv(output_folder+'/result_sustain.csv',index=False)
