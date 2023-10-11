# RLCP 
The project aims at developing a locally-weighted conformal prediction algorithm, which can attain meaningful local coverage guarantees in a model-free manner.

## Overview
We provide the exact code for reproducing the numerical results in the paper " "

## Folders
- `utils/`: contains the helper functions for all the experiments.
- `simulations/`: contains the codes for simulation experiments in the paper.
- `real_data/`: contains the codes for real data experiments in the paper.
- `results/`: contains the results from the real data and simulation experiments.

## Guide for the codes in `simulations/` folder
### Codes for reproducing results in Section 5.1
- marginal coverage `(Figure 1)`: `mc_cov.R/`
- average prediction interval `(Figure 2)`: `avg_interv.R/`
- marginal coverage on subsets with constant bandwidth `(Figure 3)`: `mc_onA_const_h.R/`
- marginal coverage on subsets with dimension-adaptive bandwidth `(Figure 4)`: `mc_onA_opt_h.R/`

### Other simulation results from Section 6 and appendix
- effect of randomization on RLCP prediction interval in simulation settings `(Figure 7 & 8 left)`: `simul_rand_effect.R/`
- coverage of mRLCP prediction interval `(Figure 9)`: `mRLCP.R/`


## Guide for the codes in `real_data/` folder
### Codes for reproducing results in Section 5.2
- marginal coverage and conditional coverage `(Figure 5 & 6)`: `real_exp.R/`
- quantitative analysis on effect of randomization on RLCP prediction interval in real data experiments `(Figure 8 right)`: `real_dev.R/`

## Real Data
For our real data experiments, we have used the `Train_Data.csv` of medical insurance dataset from [kaggle datasets](https://www.kaggle.com/datasets/rajgupta2019/medical-insurance-dataset). For reproducing the real data experiments, please download the data from this above link, and load it into your R before running the codes from `real_data/` folder.