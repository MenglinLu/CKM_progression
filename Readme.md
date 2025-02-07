# Characterizing Subtypes, Progressive Trajectories and Lifestyle Effects in Early-Stage CKM Syndrome

## Overview

This library is designed to identify subtypes of early CKM syndrome and track the progression of the disease. It also analyzes the association between the inferred subtypes and different outcomes, and proposes a model that integrates subtype and progression information into CVD prediction.  
In the inference of subtypes and stages, participants' kidney and metabolic indicators are used as input data, with the corresponding subtype and stage for each patient being the output. The code for the SuStaIn model is referenced from [https://github.com/ucl-pond/pySuStaIn](https://github.com/ucl-pond/pySuStaIn).  
In the prediction based on the constructed model, the inferred subtypes and stages, along with the indicators included in the traditional PCE model, are used as input data, with the 10-year CVD risk as the output.

## OS requirements

This package was trained and tested under Ubuntu 18.04. Modifications may be necessary to run it on other platforms.

## Python Dependencies

- `numpy`
- `pandas`
- `scikit-learn`
- `pySuStaIn`
- `Pillow`
- `tqdm`

## R Dependencies

- `ggplot2`
- `ggiraphExtra`
- `tidyverse`
- `ggalluvial`
- `survival`
- `survminer`
- `KMunicate`
- `dplyr`

## Code structure

### sustain
- `sustain.py`: Used for subtype and stage inference of early-stage CKM syndrome
- `sustain-cv.py`: Used for the selection of the optimal number of subtypes
- `sample_data.csv`: Sample dataset for input
- `simfuncs.py`: For implementation of SuStaIn model

### clinical profiles
- `biomarker statistics.py`: Used for calculating the mean values of biomarkers across different subtypes and phases
- `biomarker changes plotting.R`: Used for plotting radar charts of biomarkers by subtypes and phases
- `disease proportion statistics.py`: Used for calculating the proportion of various conditions across different subtypes and phases
- `disease progression trajectories plotting.R`: Used for plotting percentage-stacked area charts of diseases by subtypes and phases

### longitudinal validation
- `biomarker changes validation.py`: Used for analyzing observed biomarker distributions at follow-up for longitudinal validation
- `disease progression validation.py`: Used for calculating disease proportion at follow-up for longitudinal validation

### survival analysis
- `HR_subtypes.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype
- `HR_phases.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype and each phase
- `survival curves.R`: Used for visualization

### proteomic and NMR analysis
- `01_DE_prot.R`: Used to investigate the difference in proteomic profiles between subtypes and healthy controls
- `02_DE_metabo.R`: Used to investigate the difference in metabolomic profiles between subtypes and healthy controls 

### predictive analysis
- `model.py`: A 10-year CVD risk prediction model based on subtype and progression information

### lifestyle intervention
- `HR_subtype-Lifestyle.R`: Used to calculate the HR values of healthy lifestyle behaviours for the long-term risk of CVD in comparison to unhealthy counterparts for each subtype
- `HR_subtype-Lifestyle.R`: Used to calculate the HR values of healthy lifestyle behaviours for the long-term risk of CVD in comparison to unhealthy lifestyle behaviours across different subtypes and phases

## How to install and run SuStaIn model using sample data
# Characterizing Subtypes, Progressive Trajectories and Lifestyle Effects in Early-Stage CKM Syndrome

## Overview

This library is designed to identify subtypes of early CKM syndrome and track the progression of the disease. It also analyzes the association between the inferred subtypes and different outcomes, and proposes a model that integrates subtype and progression information into CVD prediction.  
In the inference of subtypes and stages, participants' kidney and metabolic indicators are used as input data, with the corresponding subtype and stage for each patient being the output. The code for the SuStaIn model is referenced from [https://github.com/ucl-pond/pySuStaIn](https://github.com/ucl-pond/pySuStaIn).  
In the prediction based on the constructed model, the inferred subtypes and stages, along with the indicators included in the traditional PCE model, are used as input data, with the 10-year CVD risk as the output.

## OS requirements

This package was trained and tested under Ubuntu 18.04. Modifications may be necessary to run it on other platforms.

## Python Dependencies

- `python`==3.8
- `numpy`
- `pandas`
- `scikit-learn`
- `pySuStaIn`
- `Pillow`
- `tqdm`

## R Dependencies

- `ggplot2`
- `ggiraphExtra`
- `tidyverse`
- `ggalluvial`
- `survival`
- `survminer`
- `KMunicate`
- `dplyr`

## Code structure

### sustain
- `sustain.py`: Used for subtype and stage inference of early-stage CKM syndrome
- `sustain-cv.py`: Used for the selection of the optimal number of subtypes
- `sample_data.csv`: Sample dataset for input
- `simfuncs.py`: For implementation of SuStaIn model

### clinical profiles
- `biomarker statistics.py`: Used for calculating the mean values of biomarkers across different subtypes and phases
- `biomarker changes plotting.R`: Used for plotting radar charts of biomarkers by subtypes and phases
- `disease proportion statistics.py`: Used for calculating the proportion of various conditions across different subtypes and phases
- `disease progression trajectories plotting.R`: Used for plotting percentage-stacked area charts of diseases by subtypes and phases

### longitudinal validation
- `biomarker changes validation.py`: Used for analyzing observed biomarker distributions at follow-up for longitudinal validation
- `disease progression validation.py`: Used for calculating disease proportion at follow-up for longitudinal validation

### survival analysis
- `HR_subtypes.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype
- `HR_phases.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype and each phase
- `survival curves.R`: Used for visualization

### proteomic and NMR analysis
- `01_DE_prot.R`: Used to investigate the difference in proteomic profiles between subtypes and healthy controls
- `02_DE_metabo.R`: Used to investigate the difference in metabolomic profiles between subtypes and healthy controls 

### predictive analysis
- `model.py`: A 10-year CVD risk prediction model based on subtype and progression information

### lifestyle intervention
- `HR_subtype-Lifestyle.R`: Used to calculate the HR values of healthy lifestyle behaviours for the long-term risk of CVD in comparison to unhealthy counterparts for each subtype
- `HR_subtype-Lifestyle.R`: Used to calculate the HR values of healthy lifestyle behaviours for the long-term risk of CVD in comparison to unhealthy lifestyle behaviours across different subtypes and phases

## How to install and run SuStaIn model using sample data
1. Create a virtual environment‌: conda create -n sustain_test python=3.8
2. Activate the virtual environment: conda activate sustain_test
3. Install the pySuStaIn package referring to https://github.com/ucl-pond/pySuStaIn: pip install git+https://github.com/ucl-pond/pySuStaIn
4. Git clone this repo: git clone https://github.com/MenglinLu/CKM_progression.git
5. Navigate to the main directory and run `pip install -r requirements.txt`
6. Navigate to the sustain directory and run `python sustain.py` to conduct subtype and stage inference using sample data
7. The inference results are saved in the sustain folder
