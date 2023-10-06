# RLCP 
The project aims at developing a locally-weighted conformal prediction algorithm, which can attain meaningful local coverage guarantees in a model-free manner.

## Overview
We provide the exact code for reproducing the numerical results in the paper " "

## Folders
- 'utils/': contains the helper functions for all the experiments.
- 'simulations/': contains the codes for simulation experiments in the paper.
- 'real_data/': contains the codes for real data experiments in the paper.
- 'results/': contains the results from the real data and simulation experiments.

## Guide for the codes in 'simulations/' folder
### Codes for reproducing results in Section 5.1
- 'marginal coverage (Figure 1)': 'mc_cov.R/'
- 'average prediction interval (Figure 2)': 'cc_coff.R/'
- 'marginal coverage on subsets with constant bandwidth (Figure 3)': 'mc_onA_const_h.R/'
- 'marginal coverage on subsets with dimension-adaptive bandwidth (Figure 4)': 'mc_onA_opt_h.R/'

### Other simulation results from Section 6 and appendix
- 'visualizing effect of randomization on RLCP prediction interval (Figure 7)': '/'
- 'quantitative analysis on effect of randomization on RLCP prediction interval in simulation settings (Figure 8 left)': '/'
- 'coverage of mRLCP prediction interval (Figure 9)': 'mRLCP.R/'


## Guide for the codes in 'real_data/' folder
### Codes for reproducing results in Section 5.2
- 'marginal coverage and conditional coverage (Figure 5 & 6)': 'real_exp.R/'
- 'quantitative analysis on effect of randomization on RLCP prediction interval in real data experiments (Figure 8 right)': '/'