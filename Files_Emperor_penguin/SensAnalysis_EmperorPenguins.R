################################################################

# Original author: Stephanie Jenouvrier
# Study DOI: 10.1111/j.1365-2486.2012.02744.x
# This code was converted from MATLAB to R
# Edited by: Esin Ickin
# Date: Nov 27 2024

#################################################################

# main code: main_proj_Stef_all_wvec.m
# don't need nSimE which are extreme events
# but need nsim parameter uncertainty
# line 51 extreme events, don't need the lines before
# nc = 66 (66 colonies) doesn't really matter
# WEC = N = pop vector
# what we need is from line 89
# MPM with 5 stages

# beta.m are the beta parameters
# covbeta.m = covariance-variance matrix
# betarand.m = function where you sample parameters from mvdist
# simsamp.m = function to account for initial conditions --> want SSD
# project pop and then get asymptotic lambda
# theta is what we need
# parameter_rand.m??

# 4 sea ice covariates in matrix with year and season
# but only strongest effects in laying and rearing season
# ignore the frequency (like the colonization factor of the trees)
# paper equation 9 & 10 --> breederav2.m (frequency parameter)

# popmat.m = pop matrix, theta (vital rates)

# line 128 in main_proj
# save the final pop model after running for 1000 years to get SSD
# wvec = lambda

# project pop with perturbations of laying and rearing sea ice covariates


# 0) Prepare session ###########################################

rm(list=ls())

library(readxl)
library(R.matlab)
library(lqmm)
library(corpcor)

setwd("/Users/maria/Dropbox/teaching/esin/emperor penguin")

# 1) Load data ##########################################

# load sea ice covariates
SIC_old=read_excel("SIC.xlsx",sheet = 1)

SIC.max.LAY=data.frame(NB=rep(0,nrow(SIC_old)),
               LAY=rep(max(SIC_old$LAY),nrow(SIC_old)),
               IN=rep(0,nrow(SIC_old)),
               REA=rep(0,nrow(SIC_old)))

SIC.min.LAY=data.frame(NB=rep(0,nrow(SIC_old)),
                   LAY=rep(min(SIC_old$LAY),nrow(SIC_old)),
                   IN=rep(0,nrow(SIC_old)),
                   REA=rep(0,nrow(SIC_old)))
  
SIC.max.REA=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(0,nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(max(SIC_old$REA),nrow(SIC_old)))

SIC.min.REA=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(0,nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(min(SIC_old$REA),nrow(SIC_old)))

max.lay=max(SIC_old$LAY)
SIC.max.LAY.cov=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(max(SIC_old$LAY),nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(SIC_old$REA[SIC_old$LAY==max.lay],nrow(SIC_old)))

min.lay=min(SIC_old$LAY)
SIC.min.LAY.cov=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(min(SIC_old$LAY),nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(SIC_old$REA[SIC_old$LAY==min.lay],nrow(SIC_old)))

max.rea=max(SIC_old$REA)
SIC.max.REA.cov=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(SIC_old$LAY[SIC_old$REA==max.rea],nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(max(SIC_old$REA),nrow(SIC_old)))

min.rea=min(SIC_old$REA)
SIC.min.REA.cov=data.frame(NB=rep(0,nrow(SIC_old)),
                       LAY=rep(SIC_old$LAY[SIC_old$REA==min.rea],nrow(SIC_old)),
                       IN=rep(0,nrow(SIC_old)),
                       REA=rep(min(SIC_old$REA),nrow(SIC_old)))

# 2) Covariates #########################################

# Here are the combinations you would use:

sens.out=NULL

sd.L_SIC=sd(SIC_old$LAY)
mean.L_SIC=mean(SIC_old$LAY)

sd.R_SIC=sd(SIC_old$REA)
mean.R_SIC=mean(SIC_old$REA)

# problem is that SICa[2] is nowhere in the code to be found?
# only the SICa[4], SICf, and SICm

# 3) Functions ###########################################

# beta.m

  # beta = slope, intercept of relationship with covariates
  # final model

  # Return of breeders; range of estimable parameters before 1990
  # Time varying
  # B 1:13 years
  # NB 1:17 years
  bR <- c(-2.11, -2.42, -1.91, -2.18, -3.01, -2.71, -2.76, -2.79, -3.24, -1.95, 
          -2.66, -3.08, -2.68, -4.54, -5.18, -5.90, -2.97, 1.59, 2.15, 2.66, 
          1.49, 1.43, 1.60, 1.07, 2.14, 1.35, 4.38, 2.50, 0.690, 0.750)

  # beta = weight; cov female; beta 0; beta 1; beta 2;
  # cov male; beta 0; beta 1; beta 2;
  # all models

  bS <- matrix(c(
    0.1380, 2, 2.70, 2.98, 0.00, 2, 2.55, 2.74, 0.00,
    0.0933, 4, 3.45, 2.28, -17.99, 4, 3.10, 0.98, -15.00,
    0.0768, 4, 3.21, 2.48, -13.32, 2, 2.55, 2.05, 0.00,
    0.0739, 2, 2.67, 2.93, 0.00, 2, 2.60, 1.26, -5.69,
    0.0622, 4, 3.01, 2.08, -10.43, 1, 2.33, 0.00, 0.00,
    0.0480, 4, 3.23, 2.25, -15.21, 2, 2.67, 0.17, -8.56,
    0.0454, 4, 2.69, 2.82, 0.00, 1, 2.34, 0.00, 0.00,
    0.0404, 4, 2.97, 4.41, 0.00, 2, 2.55, 1.72, 0.00,
    0.0319, 4, 2.81, 3.54, 0.00, 3, 2.47, 1.49, 0.00,
    0.0236, 4, 3.01, 2.07, -10.71, 4, 2.31, -0.18, 0.00,
    0.0225, 1, 2.36, 0.00, 0.00, 1, 2.25, 0.00, 0.00,
    0.0216, 4, 3.57, 7.32, 0.00, 4, 3.12, 1.59, -12.47,
    0.0213, 3, 2.53, 1.72, 0.00, 1, 2.26, 0.00, 0.00,
    0.0193, 4, 2.69, 2.92, 0.00, 1, 2.37, 0.43, 0.00,
    0.0187, 4, 3.02, 2.27, -10.25, 3, 2.43, 1.26, 0.00,
    0.0176, 4, 2.79, 3.39, 0.00, 4, 2.39, 0.46, 0.00,
    0.0175, 2, 2.59, 2.11, 0.00, 4, 2.68, 0.46, -8.77,
    0.0174, 2, 2.51, 1.68, 0.00, 4, 2.24, -0.48, 0.00,
    0.0173, 4, 2.96, 4.42, 0.00, 2, 2.61, 1.00, -3.53,
    0.0173, 4, 2.96, 4.42, 0.00, 2, 2.61, 1.00, -3.53,
    0.0313, 4, 3.00, 2.08, -10.35, 1, 2.35, 0.42, 0.00,
    0.0148, 1, 2.37, 0.00, 0.00, 2, 2.35, 1.07, 0.00,
    0.0144, 3, 2.47, 1.40, 0.00, 2, 2.35, 1.31, 0.00,
    0.0137, 2, 2.67, 2.87, -0.08, 2, 2.61, 1.22, -5.70,
    0.0125, 2, 2.54, 1.88, 0.00, 3, 2.35, 1.17, 0.00,
    0.0114, 1, 2.36, 0.00, 0.00, 1, 2.27, 0.32, 0.00,
    0.0100, 2, 2.45, 1.23, 0.00, 1, 2.24, 0.00, 0.00,
    0.0100, 1, 2.36, 0.00, 0.00, 4, 2.23, -0.59, 0.00,
    0.0089, 2, 2.52, 1.68, 0.00, 1, 2.28, 0.26, 0.00,
    0.0088, 2, 2.59, 2.11, 0.00, 4, 2.68, 0.46, -8.78,
    0.0078, 3, 2.46, 1.35, 0.00, 2, 2.45, 0.36, -4.90,
    0.0073, 1, 2.35, 0.00, 0.00, 2, 2.45, 0.07, -5.45,
    0.0068, 4, 2.38, 0.38, 0.00, 1, 2.25, 0.00, 0.00,
    0.0065, 2, 2.61, 1.89, -1.05, 4, 2.68, 0.46, -8.99,
    0.0059, 1, 2.37, 0.00, 0.00, 3, 2.32, 0.79, 0.00,
    0.0058, 3, 2.47, 1.22, 0.00, 4, 2.54, 0.27, -6.45
  ), byrow = TRUE, ncol = 9)
  
# covbeta
  covbeta <- function(nsim) {
    # BS Breeding success
    # PB Proportion of 1st return - 1st year survival
    # R Return
    # S Survival
    # W 
    # SICf
    # SICm
    # array
    
    # Breeding success ***********************************
    betaBS <- c(-0.00586821143510701, -1.75351942135837)
    covBS <- matrix(c(
      0.0484976527654650, 0,
      0, 1.71486908453073
    ), nrow = 2, byrow = TRUE)
    
    # Proportion of 1st return - 1st year survival ***********************************
    betaPB <- c(1.216946361, -2.82227397)
    covPB <- matrix(c(
      4.82142783, -1.22625655,
      -1.22625655, 0.319922
    ), nrow = 2, byrow = TRUE)
    
    # Return probabilities ***********************************
    covRmod <- matrix(c(
      -2.11,	0.23,	0.05,	0.00,	0.01,	0.01,	0.01,	0.01,	0.02,	0.01	,0.01,	0.02,	0.02,	0.01,	0.02,	0.01,	0.02,	0.01,	0.02,	0.00,	0.01,	0.01,	0.00,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -2.42,	0.18,	0.00,	0.03,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.01,	0.01,	0.01,	0.02,	0.01,	0.01,	0.00,	0.01,	0.01,	0.00,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -1.91,	0.15,	0.01,	0.01,	0.02,	0.01,	0.01,	0.01,	0.02,	0.01,	0.01,	0.02,	0.02,	0.01,	0.02,	0.01,	0.02,	0.01,	0.02,	0.00,	0.01,	0.01,	0.00,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -2.18,	0.15,	0.01,	0.01,	0.01,	0.02,	0.01,	0.01,	0.02,	0.01,	0.01,	0.02,	0.02,	0.01,	0.02,	0.02,	0.02,	0.01,	0.02,	0.00,	0.01,	0.01,	0.01,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -3.01,	0.20,	0.01,	0.01,	0.01,	0.01,	0.04,	0.01,	0.02,	0.01,	0.01,	0.02,	0.02,	0.01,	0.02,	0.02,	0.02,	0.01,	0.02,	0.00,	0.01,	0.02,	0.00,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -2.71,	0.18,	0.01,	0.01,	0.01,	0.01,	0.01,	0.03,	0.02,	0.01,	0.02,	0.02,	0.03,	0.02,	0.02,	0.02,	0.02,	0.02,	0.02,	0.00,	0.01,	0.00,	0.01,	0.01,	-0.02,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -2.76,	0.23,	0.02,	0.01,	0.02,	0.02,	0.02,	0.02,	0.05,	0.01,	0.02,	0.03,	0.04,	0.02,	0.03,	0.02,	0.03,	0.02,	0.02,	0.00,	0.01,	0.00,	0.00,	0.02,	-0.03,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
      -2.79,	0.22,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.05,	-0.01,	0.02,	0.02,	0.02,	0.02,	0.02,	0.02,	0.01,	0.02,	0.00,	0.01,	0.01,	0.01,	0.00,	0.01,	-0.01,	0.00,	0.00,	0.01,	0.01,	0.00,	0.00,
      -3.24,	0.32,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.02,	-0.01,	0.10,	0.02,	0.03,	0.02,	0.03,	0.02,	0.03,	0.02,	0.03,	0.00,	0.01,	0.01,	0.01	,0.01,	-0.02,	0.00	,0.00,	0.00,	0.01,	0.01,	0.00,	0.00,
      -1.95,	0.22,	0.02,	0.01,	0.02,	0.02,	0.02,	0.02,	0.03,	0.02,	0.02,	0.05	,0.03	,0.02,	0.03,	0.03,	0.03,	0.02	,0.03	,0.01,	0.01,	0.01,	0.01,	0.01,	-0.02,	-0.01	,0.00	,0.00	,0.01,	0.01	,0.00,	0.00,
      -2.66	,0.38,	0.02,	0.02,	0.02,	0.02	,0.02,	0.03,	0.04,	0.02,	0.03,	0.03,	0.15	,0.01,	0.04,	0.04	,0.05,	0.03,	0.04	,0.01	,0.02	,0.02	,0.01	,0.02,	-0.05,	-0.01,	0.02,	-0.01,	0.00,	0.01,	0.00,	0.00,
      -3.08,	0.29	,0.01,	0.01,	0.01	,0.01	,0.01	,0.02,	0.02,	0.02,	0.02,	0.02	,0.01,	0.08,	0.03,	0.03,	0.03,	0.03,	0.03,	0.00,	0.01,	0.01,	0.01,	0.01,	-0.01,	0.00,	0.00,	-0.01,	0.02,	0.01	,0.00	,0.00,
      -2.68,	0.24,	0.02,	0.01	,0.02	,0.02	,0.02,	0.02	,0.03,	0.02	,0.03,	0.03	,0.04,	0.03,	0.06,	0.03	,0.04	,0.03,	0.04	,0.01,	0.01,	0.02,	0.01,	0.01,	-0.02,	0.00,	0.00,	0.00,	0.02,	0.02,	0.00,	0.00,
      2.50	,0.38,	0.00,	0.00,	0.00,	0.00,	0.00	,0.00,	0.00,	0.01,	0.01,	0.01,	0.01	,0.01,	0.02,	0.02,	-0.03,	0.06,	0.04,	0.00,	0.01	,0.01,	0.01,	0.00,	0.01,	0.00,	-0.01,	0.03,	0.06,	0.15	,0.01,	0.01,
      0.69,	0.12,	0.00,	0.00,	0.00,	0.00	,0.00,	0.00,	0.00,	0.00	,0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.00,	0.00,	0.01,	0.00	,0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.02,	0.01,	0.02	,0.00,
      0.75,	0.18,	0.00,	0.00	,0.00	,0.00,	0.00,	0.00	,0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00	,0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.02,	0.01,	0.00,	0.03,
-4.54,	0.43,	0.01	,0.01,	0.01,	0.02,	0.02	,0.02	,0.02,	0.02,	0.02,	0.03,	0.04	,0.03,	0.03,	0.19,	-0.01,	0.03,	0.04,	0.00,	0.01,	0.01	,0.01	,0.01,	-0.01,	0.00,	0.00,	0.01,	0.05,	0.02	,0.01,	0.00,
-5.18,	0.87,	0.02,	0.02,	0.02,	0.02,	0.02,	0.02,	0.03,	0.02,	0.03,	0.03,	0.05,	0.03,	0.04,	-0.01,	0.75,	-0.30,	0.04,	0.01,	0.01,	0.01,	0.01,	0.01,	-0.02,	0.00	,0.01,	0.00,	0.04,	-0.03,	0.00,	0.00,
-5.90,	1.08,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.02	,0.01	,0.02,	0.02,	0.03,	0.03,	0.03	,0.03,	-0.30,	1.17,	0.04,	0.00,	0.01,	0.02,	0.01,	0.01,	-0.01,	0.00,	0.00,	0.01,	0.05,	0.06,	0.00,	0.00,
-2.97,	0.32,	0.02,	0.01,	0.02,	0.02	,0.02,	0.02,	0.02,	0.02,	0.03,	0.03,	0.04,	0.03,	0.04,	0.04,	0.04,	0.04,	0.10,	0.01,	0.02,	0.02,	0.01,	0.01,	-0.01,	0.00	,0.00,	0.02,	0.07,	0.04,	0.01,	0.00,
1.59,	0.21,	0.00,	0.00	,0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.01,	0.00,	0.01,	0.00,	0.01,	0.00,	0.01,	0.05,	0.00,	0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,
2.15,	0.30,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.00,	0.09,	0.02,	0.01,	0.01,	-0.01,	0.00,	0.00,	0.00,	0.01,	0.01,	0.00	,0.00,
2.66,	0.58,	0.01,	0.01,	0.01,	0.01,	0.02,	0.00,	0.00,	0.01,	0.01,	0.01,	0.02,	0.01,	0.02,	0.01,	0.01,	0.02,	0.02,	0.01,	0.02,	0.34,	-0.01,	0.01,	-0.01,	0.00,	0.00,	0.01	,0.02,	0.01,	0.00,	0.00,
1.49,	0.23,	0.00,	0.00,	0.00,	0.01,	0.00,	0.01,	0.00,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.00,	0.01,	-0.01,	0.05,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.01,	0.00,	0.00,
1.43,	0.25,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.02,	0.00,	0.01,	0.01,	0.02,	0.01,	0.01,	0.01,	0.01,	0.01,	0.01,	0.00,	0.01,	0.01,	0.00,	0.06,	-0.04,	0.00,	0.00,	0.00,	0.01,	0.00,	0.00,	0.00,
1.60,	0.43,	-0.01,	-0.01,	-0.01,	-0.01,	-0.01,	-0.02,	-0.03,	0.01,	-0.02,	-0.02,	-0.05,	-0.01,	-0.02,	-0.01,	-0.02,	-0.01,	-0.01,	0.00,	-0.01,	-0.01,	0.00,	-0.04,	0.18,	-0.02,	0.00,	0.01,	0.03,	0.01,	0.00,	0.00,
1.07,	0.25,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	-0.01,	-0.01,	0.00,	-0.01,	-0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00	,0.00,	-0.02,	0.06,	0.01,	0.00,	0.00,	0.00,	0.00,	0.00,
2.14,	0.53,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.02,	0.00,	0.00,	0.00,	0.01,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.28,	-0.08,	-0.03,	-0.01,	0.00,	0.00,
1.35,	0.29,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	-0.01,	-0.01,	0.00,	0.01,	0.00,	0.01,	0.02,	0.00,	0.00,	0.01,	0.00,	0.00,	0.01,	0.00,	-0.08,	0.08,	0.03,	0.03,	0.01,	0.01,
4.38,	1.28,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.01,	0.01,	0.01,	0.00,	0.02,	0.02,	0.05,	0.04,	0.05,	0.07,	0.00,	0.01,	0.02,	0.01,	0.01,	0.03,	0.00,	-0.03,	0.03,	1.64,	0.06,	0.02,	0.02), ncol = 32, byrow = TRUE)
    
    betaR=covRmod[,1]
    covR=covRmod[,3:ncol(covRmod)]
    
    covR <- make.positive.definite(covR)
    
    # Survival probabilities from the model selection
    # 36 models (rows)
    BETAS <- matrix(c(
      0.1380, 2, 2.70, 2.98, 0.00, 2, 2.55, 2.74, 0.00,
      0.0933, 4, 3.45, 2.28, -17.99, 4, 3.10, 0.98, -15.00,
      0.0768, 4, 3.21, 2.48, -13.32, 2, 2.55, 2.05, 0.00,
      0.0739, 2, 2.67, 2.93, 0.00, 2, 2.60, 1.26, -5.69,
      0.0622, 4, 3.01, 2.08, -10.43, 1, 2.33, 0.00, 0.00,
      0.0480, 4, 3.23, 2.25, -15.21, 2, 2.67, 0.17, -8.56,
      0.0454, 4, 2.69, 2.82, 0.00, 1, 2.34, 0.00, 0.00,
      0.0404, 4, 2.97, 4.41, 0.00, 2, 2.55, 1.72, 0.00,
      0.0319, 4, 2.81, 3.54, 0.00, 3, 2.47, 1.49, 0.00,
      0.0236, 4, 3.01, 2.07, -10.71, 4, 2.31, -0.18, 0.00,
      0.0225, 1, 2.36, 0.00, 0.00, 1, 2.25, 0.00, 0.00,
      0.0216, 4, 3.57, 7.32, 0.00, 4, 3.12, 1.59, -12.47,
      0.0213, 3, 2.53, 1.72, 0.00, 1, 2.26, 0.00, 0.00,
      0.0193, 4, 2.69, 2.92, 0.00, 1, 2.37, 0.43, 0.00,
      0.0187, 4, 3.02, 2.27, -10.25, 3, 2.43, 1.26, 0.00,
      0.0176, 4, 2.79, 3.39, 0.00, 4, 2.39, 0.46, 0.00,
      0.0175, 2, 2.59, 2.11, 0.00, 4, 2.68, 0.46, -8.77,
      0.0174, 2, 2.51, 1.68, 0.00, 4, 2.24, -0.48, 0.00,
      0.0173, 4, 2.96, 4.42, 0.00, 2, 2.61, 1.00, -3.53,
      0.0173, 4, 2.96, 4.42, 0.00, 2, 2.61, 1.00, -3.53,
      0.0313, 4, 3.00, 2.08, -10.35, 1, 2.35, 0.42, 0.00,
      0.0148, 1, 2.37, 0.00, 0.00, 2, 2.35, 1.07, 0.00,
      0.0144, 3, 2.47, 1.40, 0.00, 2, 2.35, 1.31, 0.00,
      0.0137, 2, 2.67, 2.87, -0.08, 2, 2.61, 1.22, -5.70,
      0.0125, 2, 2.54, 1.88, 0.00, 3, 2.35, 1.17, 0.00,
      0.0114, 1, 2.36, 0.00, 0.00, 1, 2.27, 0.32, 0.00,
      0.0100, 2, 2.45, 1.23, 0.00, 1, 2.24, 0.00, 0.00,
      0.0100, 1, 2.36, 0.00, 0.00, 4, 2.23, -0.59, 0.00,
      0.0089, 2, 2.52, 1.68, 0.00, 1, 2.28, 0.26, 0.00,
      0.0088, 2, 2.59, 2.11, 0.00, 4, 2.68, 0.46, -8.78,
      0.0078, 3, 2.46, 1.35, 0.00, 2, 2.45, 0.36, -4.90,
      0.0073, 1, 2.35, 0.00, 0.00, 2, 2.45, 0.07, -5.45,
      0.0068, 4, 2.38, 0.38, 0.00, 1, 2.25, 0.00, 0.00,
      0.0065, 2, 2.61, 1.89, -1.05, 4, 2.68, 0.46, -8.99,
      0.0059, 1, 2.37, 0.00, 0.00, 3, 2.32, 0.79, 0.00,
      0.0058, 3, 2.47, 1.22, 0.00, 4, 2.54, 0.27, -6.45
    ), ncol = 9, byrow = TRUE)
    
    # here the number of simulatio is used to weigh t
    W=BETAS[,1] # model weights
    betaS=BETAS[,c(3:5, 7:9)] # beta estimate of the regression coeff
    SICf=BETAS[,2]
    SICm=BETAS[,6] # both SIC covariates for survival female and male
    
    # % To quantify uncertainties resulting from model selection
    # % and estimation error, we used the parametric bootstrap
    # % procedure introduced by Regehr et al. (2010) and Hunter et al.
    # % (2010). To generate a bootstrap sample of a model output in
    # % this procedure, a CMR model is first selected with probability
    # % proportional to its AIC weight. 
    # %HERE arraysim describes the number of time
    # % each model is selected to draw vector of parameter
    # % values from a multivariate normal distribution with
    # % a mean equal to the estimated parameter vector and a
    # % covariance matrix obtained from the Hessian matrix of the
    # % CMR model. These covariance matrix are described by covS for each model
    # % (36 in dimension 3)
    # %
    
    arraysim <- cumsum(ceiling(nsim * W))
    arraysim[arraysim > nsim] <- 0
    array <- arraysim[arraysim > 0]
    if (tail(array, 1) != nsim) {
      array <- c(array, nsim)
    }
    
    
    # Create covS as a 3D array
    covS <- array(0, dim = c(6, 6, 36))
    
    # Populate each slice of the array
    covS[,,1] <- matrix(c(
      0.0730, 0.2820, 0, 0.0460, 0.1750, 0,
      0.2820, 1.5580, 0, 0.1670, 0.8900, 0,
      0, 0, 0, 0, 0, 0,
      0.0460, 0.1670, 0, 0.0770, 0.3260, 0,
      0.1750, 0.8900, 0, 0.3260, 1.9050, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,2] <- matrix(c(
      0.0920, 0.1540, -0.9240, -0.0360, -0.0730, 0.9450,
      0.1540, 0.9220, -0.7980, -0.0110, 0.1540, 0.9020,
      -0.9240, -0.7980, 15.5870, 0.9120, 1.7910, -15.6450,
      -0.0360, -0.0110, 0.9120, 0.0770, 0.1250, -0.9060,
      -0.0730, 0.1540, 1.7910, 0.1250, 0.8750, -1.7150,
      0.9450, 0.9020, -15.6450, -0.9060, -1.7150, 15.7260
    ), nrow = 6, byrow = TRUE)
    
    covS[,,3] <- matrix(c(
      0.2550, 0.3640, -3.3700, 0.0830, 0.3030, 0,
      0.3640, 1.3720, -3.1730, 0.1290, 0.4780, 0,
      -3.3700, -3.1730, 55.0330, -0.8890, -4.0330, 0,
      0.0830, 0.1290, -0.8890, 0.0770, 0.3190, 0,
      0.3030, 0.4780, -4.0330, 0.3190, 1.9900, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,4] <- matrix(c(
      0.0870, 0.1430, -0.2960, 0.0480, 0.1560, 0.1180,
      0.1430, 0.7240, -0.1850, 0.0720, 0.2860, 0.1700,
      -0.2960, -0.1850, 10.3420, -0.1250, -0.6050, -0.0810,
      0.0480, 0.0720, -0.1250, 0.0650, 0.0950, 0.0330,
      0.1560, 0.2860, -0.6050, 0.0950, 0.9530, 0.4000,
      0.1180, 0.1700, -0.0810, 0.0330, 0.4000, 7.6720
    ), nrow = 6, byrow = TRUE)
    
    covS[,,5] <- matrix(c(
      0.1300, 0.2180, -1.7320, 0.0510, 0.0970, 0.1510,
      0.2180, 1.1680, -1.3670, 0.0870, 0.3190, 0.1320,
      -1.7320, -1.3670, 24.5090, -0.4930, -2.0460, -0.5970,
      0.0510, 0.0870, -0.4930, 0.0650, 0.1840, 0.0150,
      0.0970, 0.3190, -2.0460, 0.1840, 1.2700, 0.6420,
      0.1510, 0.1320, -0.5970, 0.0150, 0.6420, 7.3930
    ), nrow = 6, byrow = TRUE)
    
    covS[,,6] <- matrix(c(
      0.0580, 0.1100, 0.5420, -0.0080, -0.0360, -0.5090,
      0.1100, 0.6370, 0.3970, -0.0020, -0.1310, -0.3950,
      0.5420, 0.3970, 9.8120, 0.3440, 0.6540, -9.3610,
      -0.0080, -0.0020, 0.3440, 0.0650, 0.0080, -0.3360,
      -0.0360, -0.1310, 0.6540, 0.0080, 0.8870, -0.6450,
      -0.5090, -0.3950, -9.3610, -0.3360, -0.6450, 9.6270
    ), nrow = 6, byrow = TRUE)
    
    covS[,,7] <- matrix(c(
      0.1750, 0.4180, 2.3750, 0.1100, 0.4070, -2.6790,
      0.4180, 2.0470, 2.5550, 0.2090, 0.8270, -2.8010,
      2.3750, 2.5550, 33.2860, 1.6770, 3.2020, -30.4550,
      0.1100, 0.2090, 1.6770, 0.1250, 0.2850, -1.6300,
      0.4070, 0.8270, 3.2020, 0.2850, 1.6980, -2.9270,
      -2.6790, -2.8010, -30.4550, -1.6300, -2.9270, 34.0750
    ), nrow = 6, byrow = TRUE)
    
    covS[,,8] <- matrix(c(
      0.0830, 0.1710, 0.9130, 0.0270, 0.0950, -0.7720,
      0.1710, 0.8950, 0.7460, 0.0580, 0.3150, -0.6370,
      0.9130, 0.7460, 12.7300, 0.5580, 1.4380, -11.8430,
      0.0270, 0.0580, 0.5580, 0.0650, 0.0940, -0.5260,
      0.0950, 0.3150, 1.4380, 0.0940, 1.1420, -1.2280,
      -0.7720, -0.6370, -11.8430, -0.5260, -1.2280, 12.4940
    ), nrow = 6, byrow = TRUE)
    
    covS[,,9] <- matrix(c(
      0.1470, 0.2320, 0.1060, 0.0500, 0.1690, -0.0370,
      0.2320, 1.1700, 0.3650, 0.0970, 0.4380, -0.1470,
      0.1060, 0.3650, 10.6290, 0.3540, 1.0200, -0.8530,
      0.0500, 0.0970, 0.3540, 0.0650, 0.1240, -0.3050,
      0.1690, 0.4380, 1.0200, 0.1240, 1.1450, -1.1390,
      -0.0370, -0.1470, -0.8530, -0.3050, -1.1390, 10.2140
    ), nrow = 6, byrow = TRUE)
    
    covS[,,10] <- matrix(c(
      0.1030, 0.1850, -0.9130, 0.0320, 0.1120, 0.0910,
      0.1850, 0.8510, -0.5460, 0.0710, 0.2760, 0.0930,
      -0.9130, -0.5460, 10.6820, -0.2280, -0.9450, -0.2510,
      0.0320, 0.0710, -0.2280, 0.0650, 0.1010, 0.0150,
      0.1120, 0.2760, -0.9450, 0.1010, 0.9520, 0.3410,
      0.0910, 0.0930, -0.2510, 0.0150, 0.3410, 6.8720
    ), nrow = 6, byrow = TRUE)
    
    
    covS[,,11] <- matrix(c(
      0.0910, 0.1400, 1.1540, 0.0290, 0.1070, -0.8460,
      0.1400, 0.7120, 0.4870, 0.0480, 0.2360, -0.5580,
      1.1540, 0.4870, 11.3120, 0.4510, 1.1610, -10.2370,
      0.0290, 0.0480, 0.4510, 0.0650, 0.0780, -0.4270,
      0.1070, 0.2360, 1.1610, 0.0780, 0.9110, -1.0300,
      -0.8460, -0.5580, -10.2370, -0.4270, -1.0300, 10.6590
    ), nrow = 6, byrow = TRUE)
    
    covS[,,12] <- matrix(c(
      0.1980, 0.2890, -2.1130, 0.1120, 0.4290, 0.0260,
      0.2890, 1.5180, -1.3700, 0.1770, 0.6730, 0.0540,
      -2.1130, -1.3700, 21.8020, -0.3650, -1.8720, -0.3970,
      0.1120, 0.1770, -0.3650, 0.1250, 0.2490, 0.0050,
      0.4290, 0.6730, -1.8720, 0.2490, 1.4320, 0.5310,
      0.0260, 0.0540, -0.3970, 0.0050, 0.5310, 7.4890
    ), nrow = 6, byrow = TRUE)
    
    
    covS[,,13] <- matrix(c(
      0.1250, 0.2110, 0.7350, 0.0480, 0.1640, -0.5920,
      0.2110, 1.0540, 0.5690, 0.0850, 0.3940, -0.3740,
      0.7350, 0.5690, 11.7320, 0.2430, 1.0320, -10.4620,
      0.0480, 0.0850, 0.2430, 0.0650, 0.1060, -0.2450,
      0.1640, 0.3940, 1.0320, 0.1060, 1.0410, -0.8230,
      -0.5920, -0.3740, -10.4620, -0.2450, -0.8230, 10.9140
    ), nrow = 6, byrow = TRUE)
    
    covS[,,14] <- matrix(c(
      0.0920, 0.1880, -0.7830, 0.0550, 0.1840, 0.0730,
      0.1880, 0.8840, -0.5610, 0.0720, 0.2960, 0.0850,
      -0.7830, -0.5610, 10.4510, -0.2490, -1.0950, -0.2780,
      0.0550, 0.0720, -0.2490, 0.0610, 0.1120, 0.0180,
      0.1840, 0.2960, -1.0950, 0.1120, 1.0040, 0.3910,
      0.0730, 0.0850, -0.2780, 0.0180, 0.3910, 6.8590
    ), nrow = 6, byrow = TRUE)
    
    covS[,,15] <- matrix(c(
      0.1030, 0.1670, -0.7320, 0.0670, 0.1880, 0.0850,
      0.1670, 0.7120, -0.4780, 0.0750, 0.2850, 0.0910,
      -0.7320, -0.4780, 9.9260, -0.2170, -0.9600, -0.2690,
      0.0670, 0.0750, -0.2170, 0.0520, 0.0930, 0.0200,
      0.1880, 0.2850, -0.9600, 0.0930, 0.9020, 0.3380,
      0.0850, 0.0910, -0.2690, 0.0200, 0.3380, 6.4390
    ), nrow = 6, byrow = TRUE)
    
    
    covS[,,16] <- matrix(c(
      0.1340, 0.2420, -1.0010, 0.0890, 0.2270, 0.1030,
      0.2420, 1.1470, -0.7940, 0.1220, 0.4540, 0.1360,
      -1.0010, -0.7940, 12.1570, -0.2900, -1.2980, -0.3470,
      0.0890, 0.1220, -0.2900, 0.0690, 0.1260, 0.0210,
      0.2270, 0.4540, -1.2980, 0.1260, 1.1760, 0.4190,
      0.1030, 0.1360, -0.3470, 0.0210, 0.4190, 7.1480
    ), nrow = 6, byrow = TRUE)
    
    
    covS[,,17] <- matrix(c(
      0.1210, 0.1980, -0.9250, 0.0680, 0.1920, 0.0910,
      0.1980, 0.8530, -0.5980, 0.0850, 0.3220, 0.1050,
      -0.9250, -0.5980, 11.4130, -0.2600, -1.1510, -0.3070,
      0.0680, 0.0850, -0.2600, 0.0600, 0.1120, 0.0180,
      0.1920, 0.3220, -1.1510, 0.1120, 1.0340, 0.3890,
      0.0910, 0.1050, -0.3070, 0.0180, 0.3890, 6.7390
    ), nrow = 6, byrow = TRUE)
    
    
    covS[,,18] <- matrix(c(
      0.0970, 0.1700, -0.6750, 0.0590, 0.1770, 0.0730,
      0.1700, 0.7970, -0.4970, 0.0740, 0.2810, 0.0900,
      -0.6750, -0.4970, 9.8710, -0.2150, -0.9130, -0.2470,
      0.0590, 0.0740, -0.2150, 0.0510, 0.0900, 0.0170,
      0.1770, 0.2810, -0.9130, 0.0900, 0.8390, 0.3180,
      0.0730, 0.0900, -0.2470, 0.0170, 0.3180, 6.3330
    ), nrow = 6, byrow = TRUE)
    
    covS[,,19] <- matrix(c(
      0.0860, 0.1530, -0.6180, 0.0540, 0.1590, 0.0670,
      0.1530, 0.6920, -0.4240, 0.0670, 0.2390, 0.0840,
      -0.6180, -0.4240, 8.7290, -0.1920, -0.8030, -0.2170,
      0.0540, 0.0670, -0.1920, 0.0480, 0.0850, 0.0150,
      0.1590, 0.2390, -0.8030, 0.0850, 0.7700, 0.2890,
      0.0670, 0.0840, -0.2170, 0.0150, 0.2890, 5.9870
    ), nrow = 6, byrow = TRUE)
    
    covS[,,20] <- matrix(c(
      0.1100, 0.2000, -0.8120, 0.0650, 0.1850, 0.0780,
      0.2000, 0.8940, -0.5450, 0.0800, 0.3080, 0.0950,
      -0.8120, -0.5450, 10.3520, -0.2390, -1.0320, -0.2730,
      0.0650, 0.0800, -0.2390, 0.0600, 0.1110, 0.0190,
      0.1850, 0.3080, -1.0320, 0.1110, 1.0030, 0.3750,
      0.0780, 0.0950, -0.2730, 0.0190, 0.3750, 6.4530
    ), nrow = 6, byrow = TRUE)
    
    covS[,,21] <- matrix(c(
      0.0910, 0.1740, -0.6930, 0.0590, 0.1740, 0.0740,
      0.1740, 0.8180, -0.4640, 0.0690, 0.2540, 0.0870,
      -0.6930, -0.4640, 9.4070, -0.2090, -0.8670, -0.2300,
      0.0590, 0.0690, -0.2090, 0.0490, 0.0890, 0.0160,
      0.1740, 0.2540, -0.8670, 0.0890, 0.8100, 0.3000,
      0.0740, 0.0870, -0.2300, 0.0160, 0.3000, 5.9370
    ), nrow = 6, byrow = TRUE)
    
    covS[,,22] =matrix(c(
                   0.0240,         0,         0,    0.0180,    0.0210,         0,
                   0,         0,         0,         0,         0,         0,
                   0,         0,         0,         0,         0,         0,
                   0.0180,         0,         0,    0.0390,    0.1430,         0,
                   0.0210,         0,         0,    0.1430,    1.2650,         0,
                   0,         0,         0,         0,         0,         0),nrow = 6, byrow = TRUE)

    covS[,,23] =matrix(c(
                   0.033,	0.098,	0.000,	0.019,	0.035,	0.000,
                   0.098,	0.929,	0.000,	0.019,	0.154,	0.000,
                   0.000,	0.000,	0.000,	0.000,	0.000,	0.000,
                   0.019,	0.019,	0.000,	0.040,	0.150,	0.000,
                   0.035,	0.154,	0.000,	0.150,	1.243,	0.000,
                   0.000,	0.000,	0.000,	0.000,	0.000,	0.000),nrow = 6, byrow = TRUE)

    covS[,,24] <- matrix(c(
      0.0700, 0.2540, -0.0580, 0.0370, 0.1170, -0.0320,
      0.2540, 4.7540, 11.8010, 0.1460, 1.7110, 3.7190,
      -0.0580, 11.8010, 42.8790, 0.0730, 4.0680, 14.4090,
      0.0370, 0.1460, 0.0730, 0.0650, 0.0560, -0.6270,
      0.1170, 1.7110, 4.0680, 0.0560, 3.3350, 9.1010,
      -0.0320, 3.7190, 14.4090, -0.6270, 9.1010, 39.3270
    ), nrow = 6, byrow = TRUE)
    
    covS[,,25] <- matrix(c(
      0.0460, 0.1580, 0, 0.0210, 0.0280, 0,
      0.1580, 1.0680, 0, 0.0410, 0.1870, 0,
      0, 0, 0, 0, 0, 0,
      0.0210, 0.0410, 0, 0.0360, 0.1370, 0,
      0.0280, 0.1870, 0, 0.1370, 1.3540, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,26] <- matrix(c(
      0.0240, 0, 0, 0.0160, -0.0010, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0.0160, 0, 0, 0.0240, 0.0210, 0,
      -0.0010, 0, 0, 0.0210, 0.2450, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,27] <- matrix(c(
      0.0330, 0.0970, 0, 0.0160, 0, 0,
      0.0970, 0.9660, 0, 0.0030, 0, 0,
      0, 0, 0, 0, 0, 0,
      0.0160, 0.0030, 0, 0.0220, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,28] <- matrix(c(
      0.0230, 0, 0, 0.0160, 0.0040, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0.0160, 0, 0, 0.0220, 0.0130, 0,
      0.0040, 0, 0, 0.0130, 0.8130, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,29] <- matrix(c(
      0.0440, 0.1510, 0, 0.0180, -0.0060, 0,
      0.1510, 1.0600, 0, 0.0220, -0.0410, 0,
      0, 0, 0, 0, 0, 0,
      0.0180, 0.0220, 0, 0.0240, 0.0180, 0,
      -0.0060, -0.0410, 0, 0.0180, 0.2520, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,30] <- matrix(c(
      0.0590, 0.2200, 0, 0.0510, 0.0680, -0.5740,
      0.2200, 1.3450, 0, 0.1660, 0.3000, -2.4350,
      0, 0, 0, 0, 0, 0,
      0.0510, 0.1660, 0, 0.1610, 0.2700, -2.3960,
      0.0680, 0.3000, 0, 0.2700, 1.1130, -4.4290,
      -0.5740, -2.4350, 0, -2.3960, -4.4290, 41.5970
    ), nrow = 6, byrow = TRUE)
    
    covS[,,31] <- matrix(c(
      0.033, 0.098, 0, 0.017, 0.035, 0.042,
      0.098, 0.942, 0, 0.009, 0.190, 0.321,
      0, 0, 0, 0, 0, 0,
      0.017, 0.009, 0, 0.052, 0, 0,
      0.035, 0.190, 0, 0, 1.867, 5.371,
      0.042, 0.321, 0, 0, 5.371, 29.415
    ), nrow = 6, byrow = TRUE)
    
    covS[,,32] <- matrix(c(
      0.0230, 0, 0, 0.0160, 0.0170, 0.0110,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      0.0160, 0, 0, 0.0520, -0.0380, -0.7600,
      0.0170, 0, 0, -0.0380, 1.7250, 5.0540,
      0.0110, 0, 0, -0.7600, 5.0540, 29.1130
    ), nrow = 6, byrow = TRUE)
    
    covS[,,33] <- matrix(c(
      0.0250, 0.0160, 0, 0.0150, 0, 0,
      0.0160, 0.2010, 0, -0.0020, 0, 0,
      0, 0, 0, 0, 0, 0,
      0.0150, -0.0020, 0, 0.0230, 0, 0,
      0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0
    ), nrow = 6, byrow = TRUE)
    
    covS[,,34] <- matrix(c(
      0.0610, 0.1110, -0.4170, 0.0510, 0.0650, -0.6270,
      0.1110, 2.9130, 7.6240, 0.1180, 0.3130, -0.7720,
      -0.4170, 7.6240, 34.6630, -0.1730, 0.1150, 6.8610,
      0.0510, 0.1180, -0.1730, 0.0610, 0.3700, -1.1230,
      0.0650, 0.3130, 0.1150, 0.3700, 3.2160, 10.5890,
      -0.6270, -0.7720, 6.8610, -1.1230, 10.5890, 48.3080
    ), nrow = 6, byrow = TRUE)
    
    covS[,,35] =matrix(c(
                   
                   0.0240,         0,         0,    0.0160,    0.0040,         0,
                   0,         0,         0,         0,         0,         0,
                   0,         0,        0,         0,         0,         0,                   0.0160,         0,         0,    0.0350,    0.1350,         0,
                   0.0040,         0,         0,    0.1350,    1.4900,         0,
                   0,         0   ,      0  ,       0 ,        0    ,     0), nrow = 6, byrow = TRUE)
    
    
    covS[,,36] =matrix(c(
                   
                   0.0340,   0.1030,         0,    0.0260,    0.0350,   -0.1980,
                   0.1030,    1.0160,         0,    0.0220,    0.1660,   -0.3060,
                   0,         0,         0,         0 ,        0,         0,
                   0.0260,    0.0220  ,       0 ,   0.1160  ,  0.1930,   -1.7220,
                   0.0350,    0.1660  ,       0,    0.1930,    1.0080,   -3.2560,
                   -0.1980,   -0.3060,         0,   -1.7220,   -3.2560,   31.7880),nrow = 6, byrow = TRUE)

    # Returning all variables
    list(
      betaBS=betaBS, covBS=covBS, betaPB=betaPB, covPB=covPB, betaR=betaR, covR=covR, betaS=betaS, covS=covS, W=W ,SICf=SICf, SICm=SICm, array=array
    )
  }



# betarand
betarand <- function(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS) {
  # Breeding success
  bBS <- MASS::mvrnorm(1, mu = betaBS, Sigma = covBS)
  
  # Proportion of 1st return - 1st year survival
  bPB <- MASS::mvrnorm(1, mu = betaPB, Sigma = covPB)
  
  # Return probabilities
  bR <- MASS::mvrnorm(1, mu = betaR, Sigma = covR) # Time-varying estimates
  
  # Survival probabilities
  bS <- matrix(0, nrow = nrow(betaS), ncol = 6)
  for (j in 1:nrow(betaS)) { # Weight
    bS[j, ] <- MASS::mvrnorm(1, mu = betaS[j, ], Sigma = make.positive.definite(covS[, , j]))
    if (bS[j, 3] > 0) bS[j, 3] <- 0
    if (bS[j, 6] > 0) bS[j, 6] <- 0
  }
  
  # Return as a list
  list(bBS = bBS, bPB = bPB, bR = bR, bS = bS)
}

# breederav2
breederav2 <- function(THETA, nt) {
  
  # % calculate the number of males Mav and Females Fav avalaible for mating
  # % at time t+1
  # % ASR: proportion of males in the breeders class at time t
  # %          1 2  3   4   5   6  7   8  9
  # %   THETA=[s;BS;Ppb;Pnb;Pb;Spb;Sf;Sm;S0];
  # 
  
  
  # Initialize matrices to store the number of females (Fav) and males (Mav) available for mating
  Fav <- matrix(0, nrow = 5, ncol = nt)
  Mav <- matrix(0, nrow = 5, ncol = nt)
  
  # Calculate the available females and males
  Fav[1, ] <- THETA[3, ] * THETA[5, ]
  Fav[2, ] <- THETA[6, ] * THETA[8, ]
  Fav[5, ] <- 0.5 * THETA[7, ] * THETA[8, ]
  
  Mav[3, ] <- THETA[3, ] * THETA[5, ]
  Mav[4, ] <- THETA[6, ] * THETA[9, ]
  Mav[5, ] <- 0.5 * THETA[7, ] * THETA[9, ]
  
  # Return both matrices
  return(list(Fav = Fav, Mav = Mav))
}



# parameter_RAND3
parameter_RAND3 <- function(xA, SICa, bR, bS, SICf, SICm, W, tlimit) {
  # Initialize the theta matrix
  theta <- matrix(0, nrow = 9, ncol = tlimit)
  
  # Define parameters
  theta[1, ] <- 0.5 # Sex-ratio
  theta[3, ] <- 0.0561 / (1 - xA) # Proportion of 1st return
  theta[4, ] <- 0.7715 / (1 - xA) # 1st year at sea survival
  
  # Stochasticity due to interannual variability
  theta[2, ] <- invlogit(-0.0059 - 1.7535 * SICa[1,4] + rnorm(tlimit, mean = 0, sd = 1.1646)) # BS
  theta[6, ] <- invlogit(bR[ceiling(17 * runif(tlimit))]) # Proportion of return NB
  theta[7, ] <- invlogit(bR[17 + ceiling(13 * runif(tlimit))]) # Proportion of return B
  
  # Uncomment this section if you want no environmental stochasticity
  # theta[2, ] <- invlogit(-0.0059 - 1.7535 * SICa[, 4]) # BS
  # theta[6, ] <- 0.0446 # Proportion of return NB
  # theta[7, ] <- 0.8619 # Proportion of return B
  
  # Survival uncertainties in model selection - survival is a weighted average (wAIC)
  
  # Female survival
  sf <- numeric(tlimit) # Initialize female survival vector
  for (j in seq_len(nrow(bS))) {
    sf <- sf + invlogit(
      bS[j, 1] + 
        (bS[j, 2] * SICa[1,SICf[j]]) + 
        (bS[j, 3] * SICa[1,SICf[j]]^2)
    ) * W[j]
  }
  theta[8, ] <- sf
  
  # Male survival
  sm <- numeric(tlimit) # Initialize male survival vector
  for (j in seq_len(nrow(bS))) {
    sm <- sm + invlogit(
      bS[j, 4] + 
        (bS[j, 5] * SICa[1,SICm[j]]) + 
        (bS[j, 6] * SICa[1,SICm[j]]^2)
    ) * W[j]
  }
  theta[9, ] <- sm
  
  # Pre-breeders survival (average of female and male survival)
  theta[5, ] <- (theta[9, ] + theta[8, ]) / 2
  
  return(theta)
}

# simpsamp
simpsamp <- function(k) {
  # H Caswell
  # function c = simpsamp(k)
  # returns a set of k weights uniformly sampled from the k-simplex
  # sampling by broken stick method
  
  breaks <- c(0, runif(k - 1), 1)
  breaks <- sort(breaks)
  
  intervals <- breaks[2:(k + 1)] - breaks[1:k]
  
  c <- intervals
  return(c)
}


# population matrix
popmat <- function(theta, uf, um) {
  # theta
  # 1   2   3    4   5    6    7  8   9
  # [rho, BS, Rpb, S0, Spb, Rnb, Rb, Sf, Sm]
  
  # m <- (0.5 * theta[2]) / ((theta[8]^(8 / 12)) * (theta[9]^(8 / 12)))
  m <- (theta[2]) / ((theta[8]^(8 / 12)) * (theta[9]^(8 / 12))) # probability that a breeding pair raises offspring
  
  f1 <- (1 - theta[1]) * m
  f3 <- theta[1] * m
  
  # Initialize the population matrix A with zeros
  A <- matrix(0, nrow = 5, ncol = 5)
  
  # Initialize A as a 5x5 matrix with zeros
  A <- matrix(0, nrow = 5, ncol = 5)
  
  # Assign values to the elements of A based on the provided formulas
  A[1, 1] <- (1 - uf * theta[3]) * theta[5]
  A[1, 2] <- 0
  A[1, 3] <- 0
  A[1, 4] <- 0
  A[1, 5] <- (1 - uf * theta[3]) * theta[4] * f1
  
  A[2, 1] <- 0
  A[2, 2] <- (1 - uf * theta[6]) * theta[8]
  A[2, 3] <- 0
  A[2, 4] <- 0
  A[2, 5] <- 0.5 * (1 - uf * theta[7]) * theta[8]
  
  A[3, 1] <- 0
  A[3, 2] <- 0
  A[3, 3] <- (1 - um * theta[3]) * theta[5]
  A[3, 4] <- 0
  A[3, 5] <- (1 - um * theta[3]) * theta[4] * f3
  
  A[4, 1] <- 0
  A[4, 2] <- 0
  A[4, 3] <- 0
  A[4, 4] <- (1 - um * theta[6]) * theta[9]
  A[4, 5] <- 0.5 * (1 - um * theta[7]) * theta[9]
  
  A[5, 1] <- uf * theta[3] * theta[5]
  A[5, 2] <- uf * theta[6] * theta[8]
  A[5, 3] <- um * theta[3] * theta[5]
  A[5, 4] <- um * theta[6] * theta[9]
  A[5, 5] <- 0.5 * uf * theta[7] * theta[8] + 0.5 * um * theta[7] * theta[9] + 
    uf * theta[3] * theta[4] * f1 + um * theta[3] * theta[4] * f3
  
  return(A)
}


invlogit <- function(x) {
  y <- exp(x) / (1 + exp(x))
  return(y)
}

# 4) Run simulations #################################

# SIMULATIONS

# nt = number of years = 1000
# nsimulation = number of simulations = 100
nt=1000
nsimulation=100
nsim=100
# start main code main_proj...

# Get covariance matrix of demographic parameters
  # Uncertainties on the estimated parameters
  bR <- bR # Estimate CMR, return, and survival
  bS <- bS
  xA <- 0.058  # Tag loss
  beta_results <- covbeta(nsim)
  betaBS <- beta_results$betaBS
  covBS <- beta_results$covBS
  betaPB <- beta_results$betaPB
  covPB <- beta_results$covPB
  betaR <- beta_results$betaR
  covR <- beta_results$covR
  betaS <- beta_results$betaS
  covS <- beta_results$covS
  W <- beta_results$W
  SICf <- beta_results$SICf
  SICm <- beta_results$SICm
  array <- beta_results$array
  
############# MAXIMUM LAY
  Rtot <- matrix(0, nt, nsimulation)
  WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
  x <- 1
  
  for (s in 1:nsimulation) {  # Sampling parameter uncertainty
    
    # Parameter uncertainty
    betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
    bBS <- betas$bBS
    bPB <- betas$bPB
    bR <- betas$bR
    bS <- betas$bS
    
    # POP INITIAL CONDITIONS
    r <- rep(0, nt)
    w <- array(0, dim = c(5, nt))
    
    # SAD
    wvec <- simpsamp(5)
    
    # Vital rates sampling stochasticity
    THETA <- parameter_RAND3(xA, SIC.max.LAY, bR, bS, SICf, SICm, W, nt)
    # ENS[[ens]]$SICa[, , col]
    
    # Population projection
    Fav <- breederav2(THETA, nt)$Fav
    Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
    
    for (t in 1:nt) {  # From 1900 to 2100
      # Calculate demographic parameters
      uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
      um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
      
      # Project population
      wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
      r[t] <- sum(wvec)  # Calculate growth rate
      wvec <- wvec / r[t]
      w[, t] <- wvec
    }  # Time
    
    Rtot[, s] <- r  # R for each simulation for each t for one ensemble
    WVEC[, , s] <- w
    x <- x + 1
    # Demographic sim
    
    # 5) Pop model #############
    
    # isolate the population matrix from the loop just to see if that's a problem
    
    # # nt = number of years = 1000
    # # nsimulation = number of simulations = 100
    # nt=100 #10 works but 1000 doesn't work
    # nsimulation=1
    # nsim=100
    
    # POP INITIAL CONDITIONS
    r <- rep(0, nt)
    w <- array(0, dim = c(5, nt))
    
    # SAD
    wvec <- simpsamp(5)
    
    
    THETA <- parameter_RAND3(xA, SIC.max.LAY, bR, bS, SICf, SICm, W, nt)
    
    # Population projection
    Fav <- breederav2(THETA, nt)
    Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
    
    
    
    for (t in 1:nt) {  # From 1900 to 2100
      # Calculate demographic parameters
      uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
      um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
      
      # Project population
      wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
      r[t] <- sum(wvec)  # Calculate growth rate
      wvec <- wvec / r[t]
      w[, t] <- wvec
    }  # Time
    
    Rtot[, s] <- r  # R for each simulation for each t for one ensemble
    WVEC[, , s] <- w
    x <- x + 1
  }  
  

r_fin_max_lay=Rtot[t,]


hist(r_fin_max_lay) 


############# MINIMUM LAY
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.min.LAY, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.min.LAY, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_min_lay=Rtot[t,]


hist(r_fin_min_lay) 


############# MAXIMUM REA
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.max.REA, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.max.REA, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_max_rea=Rtot[t,]


hist(r_fin_max_rea) 


############# MINIMUM REA
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.min.REA, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.min.REA, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_min_rea=Rtot[t,]


hist(r_fin_min_rea) 
        
############# MAXIMUM LAY COV
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.max.LAY.cov, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.max.LAY.cov, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_max_lay_cov=Rtot[t,]


hist(r_fin_max_lay_cov) 


############# MINIMUM LAY COV
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.min.LAY.cov, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.min.LAY.cov, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_min_lay_cov=Rtot[t,]


hist(r_fin_min_lay_cov) 


############# MAXIMUM REA
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.max.REA.cov, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.max.REA.cov, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_max_rea_cov=Rtot[t,]


hist(r_fin_max_rea_cov) 


############# MINIMUM REA
Rtot <- matrix(0, nt, nsimulation)
WVEC <- array(0, dim = c(5, nt, nsimulation))  # To save SAD
x <- 1

for (s in 1:nsimulation) {  # Sampling parameter uncertainty
  
  # Parameter uncertainty
  betas <- betarand(betaBS, covBS, betaPB, covPB, betaR, covR, betaS, covS)
  bBS <- betas$bBS
  bPB <- betas$bPB
  bR <- betas$bR
  bS <- betas$bS
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  # Vital rates sampling stochasticity
  THETA <- parameter_RAND3(xA, SIC.min.REA.cov, bR, bS, SICf, SICm, W, nt)
  # ENS[[ens]]$SICa[, , col]
  
  # Population projection
  Fav <- breederav2(THETA, nt)$Fav
  Mav <- breederav2(THETA, nt)$Mav # Calculate males and females available to mate
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Fav[, t] %*% wvec)
    um <- min(Fav[, t] %*% wvec, Mav[, t] %*% wvec) / (Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
  # Demographic sim
  
  # 5) Pop model #############
  
  # isolate the population matrix from the loop just to see if that's a problem
  
  # # nt = number of years = 1000
  # # nsimulation = number of simulations = 100
  # nt=100 #10 works but 1000 doesn't work
  # nsimulation=1
  # nsim=100
  
  # POP INITIAL CONDITIONS
  r <- rep(0, nt)
  w <- array(0, dim = c(5, nt))
  
  # SAD
  wvec <- simpsamp(5)
  
  
  THETA <- parameter_RAND3(xA, SIC.min.REA.cov, bR, bS, SICf, SICm, W, nt)
  
  # Population projection
  Fav <- breederav2(THETA, nt)
  Mav <- breederav2(THETA, nt) # Calculate males and females available to mate
  
  
  
  for (t in 1:nt) {  # From 1900 to 2100
    # Calculate demographic parameters
    uf <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Fav$Fav[, t] %*% wvec)
    um <- min(Fav$Fav[, t] %*% wvec, Mav$Mav[, t] %*% wvec) / (Mav$Mav[, t] %*% wvec)
    
    # Project population
    wvec <- popmat(THETA[, t], uf, um) %*% wvec  # Define population matrix
    r[t] <- sum(wvec)  # Calculate growth rate
    wvec <- wvec / r[t]
    w[, t] <- wvec
  }  # Time
  
  Rtot[, s] <- r  # R for each simulation for each t for one ensemble
  WVEC[, , s] <- w
  x <- x + 1
}  


r_fin_min_rea_cov=Rtot[t,]


hist(r_fin_min_rea_cov)  


sens.out=rbind(sens.out,data.frame(spec.driver="SIC laying season",
                                   driver="Temperature",
                                   driver.type="C",
                                   stage.age="all",vital.rates="all",
                                   sens= abs(r_fin_max_lay-r_fin_min_lay)/abs((max.lay-min.lay)/sd.L_SIC),
                                   l_ratio=abs(log(r_fin_max_lay/r_fin_min_lay)),
                                   cov=0,
                                   dens=0,
                                   sim=1:nsim))

sens.out=rbind(sens.out,data.frame(spec.driver="SIC laying season",
                                   driver="Temperature",
                                   driver.type="C",
                                   stage.age="all",vital.rates="all",
                                   sens= abs(r_fin_max_lay_cov-r_fin_min_lay_cov)/abs((max.lay-min.lay)/sd.L_SIC),
                                   l_ratio=abs(log(r_fin_max_lay_cov/r_fin_min_lay_cov)),
                                   cov=1,
                                   dens=0,
                                   sim=1:nsim))

sens.out=rbind(sens.out,data.frame(spec.driver="SIC rearing season",
                                   driver="Temperature",
                                   driver.type="C",
                                   stage.age="all",vital.rates="all",
                                   sens= abs(r_fin_max_rea-r_fin_min_rea)/abs((max.rea-min.rea)/sd.R_SIC),
                                   l_ratio=abs(log(r_fin_max_rea/r_fin_min_rea)),
                                   cov=0,
                                   dens=0,
                                   sim=1:nsim))

sens.out=rbind(sens.out,data.frame(spec.driver="SIC rearing season",
                                   driver="Temperature",
                                   driver.type="C",
                                   stage.age="all",vital.rates="all",
                                   sens= abs(r_fin_max_rea_cov-r_fin_min_rea_cov)/abs((max.rea-min.rea)/sd.R_SIC),
                                   l_ratio=abs(log(r_fin_max_rea_cov/r_fin_min_rea_cov)),
                                   cov=1,
                                   dens=0,
                                   sim=1:nsim))
sens_emp=data.frame(species="Aptenodytes forsteri", study.doi="10.1111/j.1365-2486.2012.02744.x",year.of.publication=2012,
                    group="Birds",continent="Antarctic",driver=sens.out$driver,driver.type="C",
                    stage.age="all",vital.rates="all",sens=sens.out$sens,cov=sens.out$cov,mat=3,n.vr=8,n.pam=42,dens=0,
                    biotic_interactions=0,lambda.sim=0,study.length=29,l_ratio=sens.out$l_ratio)


write.csv(sens_emp,"sens_emperor_peng.csv",row.names = F)

library(ggplot2)

ggplot(sens.out, aes(x=sens, color=factor(cov))) +
  geom_density()+
  facet_grid(driver~.,scales = "free")

ggplot(sens.out, aes(x=l_ratio, color=factor(cov))) +
  geom_density()+
  facet_grid(driver~.,scales = "free")

