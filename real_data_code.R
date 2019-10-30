##########################################################################
library(dplyr)
library(mice)
library(numDeriv)
library(rootSolve)
# Load the relevant code
source('semi_sensitivity_binary.R')

###############################################

# Construct a 95% CI of the average treatment effect 
# with a working model on U specified by u_vec and
# u_prob, and two sensitivity parameters delta and lambda.

obtain_mu_sd_one_dt <- function(delta, lambda, u_vec, u_prob, dataset){
  
  # Fit an outcome model and a propensity score model
  # as an initialization.
  dt_real_small = dataset
  n = dim(dt_real_small)[1]
  
  # Fit the propensity score model
  kappa_0 = glm(abd ~.-1, family = 'binomial', offset = 2*rbinom(n, 1, 0.5),  data = dt_real_small[,1:6])
  kappa = coef(kappa_0)
  
  # Fit an outcome regression model
  beta_0 = glm(vote05 ~.-1, family = 'binomial', data = dt_real_small)
  beta = coef(beta_0)
  
  # Solve the estimating equation and find the root
  root = multiroot(estimating_eq_2, c(beta,kappa), q = 5, lambda = lambda, delta = delta,
                   u_vec = u_vec, u_prob = u_prob, dataset = dt_real_small)
  param_0 = root$root
  est = param_0[6]
  
  # Compute the robust sandwich estimator of the variance
  var_beta = sandwich_variance(param_0 = param_0, q = 5, lambda = lambda, delta = delta,
                               u_vec = u_vec, u_prob = u_prob, dataset = dt_real_small)
  sd = sqrt(var_beta[6,6])
  return(c(est, sd))
}

# Do once
dataset = read.csv('dt_real_small_1.csv')
obtain_mu_sd_one_dt(0.3, 0.3, c(0,1), c(0.5, 0.5), dataset)

#########################################################
# do five datasets and pool together results
mi_mu_sd <- function(delta, lambda, u_vec, u_prob){
  cat(1, '\n')
  dataset_1 = read.csv('dt_real_small_1.csv')
  res_1 = obtain_mu_sd_one_dt(lambda, delta, u_vec, u_prob, dataset_1)  
  cat(2, '\n')
  dataset_2 = read.csv('dt_real_small_2.csv')
  res_2 = obtain_mu_sd_one_dt(lambda, delta, u_vec, u_prob, dataset_2)  
  cat(3, '\n')
  dataset_3 = read.csv('dt_real_small_3.csv')
  res_3 = obtain_mu_sd_one_dt(lambda, delta, u_vec, u_prob, dataset_3)  
  cat(4, '\n')
  dataset_4 = read.csv('dt_real_small_4.csv')
  res_4 = obtain_mu_sd_one_dt(lambda, delta, u_vec, u_prob, dataset_4)  
  cat(5, '\n')
  dataset_5 = read.csv('dt_real_small_5.csv')
  res_5 = obtain_mu_sd_one_dt(lambda, delta, u_vec, u_prob, dataset_5)  
  res = rbind(res_1, res_2, res_3, res_4, res_5)
  point_est = res[,1]
  variance_est = res[,2]^2
  mi_point_est = mean(point_est)
  total_var = mean(variance_est) + var(point_est)
  return(c(mi_point_est, sqrt(total_var)))
}

# Make Table 4 in the paper

# Binary U case
mi_mu_sd(0.2, 0.2, c(0,1), c(0.5, 0.5))
mi_mu_sd(0.5, 0.5, c(0,1), c(0.5, 0.5))
mi_mu_sd(0.6, 0.6, c(0,1), c(0.5, 0.5))
mi_mu_sd(0.7, 0.7, c(0,1), c(0.5, 0.5))
mi_mu_sd(0.8, 0.8, c(0,1), c(0.5, 0.5))

# Continuous U case
mi_mu_sd(0.2, 0.2, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(0.5, 0.5, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(0.7, 0.7, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(0.8, 0.8, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(0.9, 0.9, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(1, 1, seq(0, 1, 0.2), rep(1/6, 6))
mi_mu_sd(1.1, 1.1, seq(0, 1, 0.2), rep(1/6, 6))
