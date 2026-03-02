# Carriers-epidemic-treatment-failure
This repository contains the complete simulation code  used in the manuscript: 

Modelling the effects of treatment failure on the minor outbreak duration for carrier-related infectious disease

The code reproduces all figures and numerical results reported in the paper.
## Main Files
- centraltendency.m        # Computation of central tendecy, SD, outlier
- simdata.m                # Stochastic Simulation 
- pdfpaper2.m              # Computation of PDF for minor outbreak duration and estimate probability of extinction

## Input and Output
- Default parameter values are specified at the beginning of each script (see the values in the manuscript).
- The input parameters of simdata are beta, f, Time (final time), and P_ext_ana (analytic probability of extinction).
- The output of simdata is T_allend (simulated dataset)
- The output of centraltendency are mean, median, SD, outlier, oir (outlier influence ratio).
- The baseline critical threshold is f_c = 0.37.
- The input parameters of pdfpaper2 are all parameter values.
- The output of pdfpaper2 are PDF graph and probability of extinction when starting with carriers and symptomatic infection

## Instructions:
- Run: centraltendecy.m to reproduce Figures 6-8
- Run: pdfpaper2.m to reproduce Figures 2-5
- Simulated data can be used to perform further calculations
