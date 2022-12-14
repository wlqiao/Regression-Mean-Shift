# Space partitioning and regression maxima seeking via a mean-shift-inspired algorithm

This repository contains the R codes for the following article: 

Space partitioning and regression maxima seeking via a mean-shift-inspired algorithm.  
Qiao, W. and Shehu. A. (2022). Electronic Journal of Statistics, 16(2), 5623-5658.

## Abstract

The mean shift (MS) algorithm is a nonparametric method used to cluster sample points and find the local modes of kernel density estimates, using an idea based on iterative gradient ascent. In this paper we develop a mean-shift-inspired algorithm to estimate the maxima of regression functions and partition the sample points in the input space. We prove convergence of the sequences generated by the algorithm and derive the rates of convergence of the estimated local maxima for the underlying regression model. We also demonstrate the utility of the algorithm for data-enabled discovery through an application on biomolecular structure data.

## Description

- "RMS_funs.r": collection of functions for regression mean shift algorithm

- "RMS_simulation.r": R code for simulations (Section 4.1)

- "Simulation_result_aggregation.r": R code for Figure 4.1

- "RMS_example.r": R code for Figure 4.2

- "RMS_landscape.r": R code for Figure 4.3

- "RMS_Malaria.r": R code for Figure 4.4

- "RMS_example_NW.r": R code to show the failure of convergence if the Nadaraya-Watson regression estimator and its gradient are directly used with a mean shift idea (Remark 2.1)

## Citation

If you find this repository useful, please cite the following paper: 

> article{qiao2022space,  
  title={Space partitioning and regression maxima seeking via a mean-shift-inspired algorithm},  
  author={Qiao, Wanli and Shehu, Amarda},  
  journal={Electronic Journal of Statistics},  
  volume={16},  
  number={2},  
  pages={5623--5658},  
  year={2022},  
  publisher={Institute of Mathematical Statistics and Bernoulli Society}  
}
