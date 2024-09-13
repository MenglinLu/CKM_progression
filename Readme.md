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
- `sklearn`
- `pySuStaIn`
- `pillow`
- `argparse`
- `tqdm`

## Code structure

### sustain
- `sustain.py`: Used for subtype and stage inference of early CKM syndrome
- `sustain-cv.py`: Used for the selection of the optimal number of subtypes
- `MyDataset.csv`: Sample dataset for input

### survival analysis
- `HR_subtypes.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype
- `HR_stages.R`: Used to calculate the HR values for the long-term risk of different outcomes for each subtype and each phase
- `survival curves.R`: Used for visualization

### proteomic and NMR analysis
- `01_DE_prot.R`: Used to investigate the difference in proteomic profiles between subtypes and healthy controls
- `02_DE_metabo.R`: Used to investigate the difference in metabolomic profiles between subtypes and healthy controls 

### predictive analysis
- `model.py`: A 10-year CVD risk prediction model based on subtype and progression information
