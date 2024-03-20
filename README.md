# Data and Code for "*Constructing coarse-grained models with physics-guided Gaussian process regression*"  

## Overview 

The data and code shared here in the repository are to support the findings in the following paper titled "*Constructing coarse-grained models with physics-guided Gaussian process regression*" by Yating Fang, Qian Qian Zhao, Ryan B. Sills, and Ahmed Aziz Ezzat. It is meant to promote transparency and reproducibility in the research, rather than to be a plug-and-play software solution. Even though you need to adjust the code in order to adapt it to your own application, the concept and methodology as outlined in the paper is broadly applicable.

## Paper Abstract
Coarse-grained models describe the macroscopic mean response of a process at large scales which derives from stochastic processes at small scales. Common examples include accounting for velocity fluctuations in a turbulent fluid flow model and cloud evolution in climate models. Most existing techniques for constructing coarse-grained models feature ill-defined parameters whose values are arbitrarily chosen (e.g., a window size), are narrow in their applicability (e.g., only applicable to time series or spatial data), or cannot readily incorporate physics information. Here we introduce the concept of physics-guided Gaussian process regression as a machine-learning-based coarse-graining technique which is broadly applicable and amenable to input from known physics-based relationships. Using a pair of case studies derived from molecular dynamics simulations, we demonstrate the attractive properties and superior performance of physics-guided Gaussian processes for coarse-graining relative to prevalent benchmarks. The key advantage of Gaussian-process-based coarse-graining is its ability to seamlessly integrate data-driven and physics-based information.

## Repository Content
Case Study 1: 
* traction_data.csv (data)
* case_study_1.R (code)

Case Study 2: 
* nucleation_data.xlsx (data)
* case_study_2.R (code)

## Prerequisite
RStudio 

## Usage Instruction
Depending on the case study that you are interested in, please follow these steps to reproduce some of the results shown in paper:
1. Download the R code file and the corresponding data file
2. Ensure that the two files are in the same folder on your computer 
3. Open the code file in RStudio
4. Press the "source" button to run the entire code file
5. You may need to install some libraries if prompted by RStudio

## Code Outputs
There are mainly two outputs for either of the code files:
* A graph in the "plot window" that visualizes the training fit
* An NRMSD value printed out in the console

Interpretation of the outputs: 
* A training fit is generally considered good if it follows closely with the data, showing the corresponding model's capability to capture any trend present in the data
* A lower NRMSD is always considered good, as it indicates that the model can generalize well to unseen data

## Citation Information

If you would like to reference this repository or any of its contents, please cite as follows (to be updated upon publication): 

Fang, Y.; Zhao, Q. Q.; Sills, R. B.; Ezzat, A. A. Constructing Coarse-Grained Models with Physics-Guided Gaussian Process Regression. Manuscript under review, 2023. 

## Contact Information
Please contact Yating Fang at yf245@rutgers.edu for any questions.


