# Progress towards the Sustainable Development Goals:  a systematic analysis for estimating global, regional, and national causes of under-five mortality during 2000-2019

Code and functions for paper by Jamie PERIN, Amy MULICK, Diana YEUNG, Francisco VILLAVICENCIO, Gerard LOPEZ, Kathleen L STRONG, David PRIETO-MERINO, Simon COUSENS, Robert E BLACK and Li LIU, published at The Lancet Child & Adolescent Health, 2021.

This repository contains two folders (neonatal_cod, postneonatal_cod), each of which is organised into subfolders containing code, data and results.

## Neonatal (0-28 days) cause-of-death fractions

To reproduce the model-based estimates, the following five files located in neonatal_cod/code should be run in sequence:

1. Estimate_VA.R
2. Estimate_VR1.R
3. Estimate_VR2.R
4. Predict.R
5. Compile_VAVR.R

Functions, statistical models and datasets required to run them are included in the subfolders. Output files and datasets are also included for those wishing to view the results without running the code, except for the initial Bayesian model estimates (Estimate_XX.R, numbers 1-3 above). The files produced by these scripts are too large to store on Github.


## Child (1-59 months) cause-of-death fractions

To reproduce the model-based estimates, the following five files located in postneonatal_cod/code should be run in sequence:

1.  Estimate_VA.R
2.  Estimate_VR.R
3.  Predict_VA.R
4.  Predict_VR.R
5.  Model_Averaging.R

Functions, statistical models and datasets required to run them are included in the subfolders. Output files and datasets are also included for those wishing to view the results without running the code, except for the initial Bayesian model estimates (Estimate_XX.R, numbers 1-2 above). The files produced by these scripts are too large to store on Github.


## Combine estimated cause fractions and envelopes

The code to reproduce the final estimates is in the main file:

child_cod_2000-2019.R 

This script further parcels out and/or adjusts estimates for measles, TB, crisis, HIV, and malaria (outside high endemicity areas) deaths as described in the paper.

