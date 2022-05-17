#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_sharp.R")
source("../utils/methods_sharp_ci.R")
source("nonlinear_model.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 500         # number of observations
p = 100          # number of variables (should be even)
signal_mean = 2  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
batch <- 1
batch.size <- 1

seed = 2022
seed.model = 2022

## Parameters for P(X)
rho = 0.5
p1 = round(p/10)      # Number of treatments
p0 = p - p1           # Number of covariates

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates



#######################
## Data distribution ##
#######################

sample_X = function(n, p, p1) {
  ## Number of non-treatment variables
  p0 = p - p1
  stopifnot(p1 <= p0)
  ## Generate half of the variables from a multivariate normal distribution
  mu = rep(0,p0)
  Sigma = toeplitz(rho^(0:(p0-1)))
  X1 = matrix(rnorm(n*p0),n) %*% chol(Sigma)
  X1 = pnorm(X1-1)
  # Generate the second half of the variables from Bernoulli conditional on the first half
  U2 = matrix(runif(p1*n), n, p1)
  X2 = (U2 <= X1[,1:p1])
  # Combine the variables
  X = cbind(X2, X1)
  return(X)
}

sample_effects = function(p, p1, p_causal, p_causal_covar, signal_mean, signal_std) {
  p0 = p - p1
  nonzero_treat = sort(sample(p1, p_causal))
  nonzero_covar = sort(p1+sample(p0, p_causal_covar))
  nonzero = sort(c(nonzero_treat, nonzero_covar))
  # Flip signs for covariates
  signs = c(rep(1,p1), 2*rbinom(p0,1,0.5)-1)
  beta_center = signal_mean * (1:p %in% nonzero) * signs
  beta_center = matrix(rep(beta_center,each=n), n, p)
  beta_matrix = beta_center + matrix(rnorm(n*p, mean=0, sd=signal_std), n, p)
  beta_matrix[,-nonzero] = 0
  beta_matrix = t(beta_matrix)
  return(beta_matrix)
}

sample_Y_linear = function(X, beta, interaction_strength=0) {
  n = nrow(X)
  p = ncol(X)
  y = rep(0,n)
  for(i in 1:n) {
    y[i] = X[i,,drop=FALSE] %*% beta[,i,drop=FALSE]
  }
  if(interaction_strength>0) {
    beta_avg = rowMeans(beta)[1:p1]
    beta_support = sort(which(abs(beta_avg) > 0))
    pairs = lapply(seq(1,round(p1/2)), function(tmp) sample(beta_support, 2, replace=FALSE))
    for(pair in pairs) {
      y = y + interaction_strength * (X[,pair[1]] * X[,pair[2]] - mean( X[,pair[1]] * X[,pair[2]]) )
    }
  }
  y = y + rnorm(n)
  return(y)
}

sample_Y_tree = function(X, beta, tree, interaction_strength=0) {
  n = nrow(X)
  p = ncol(X)
  y = rep(0,n)
  df = as_tibble(X[,seq(p1+1,p)])
  y = predict(tree, df)
  for(i in 1:n) {
    y[i] = y[i] + X[i,seq(1,p1),drop=FALSE] %*% beta[seq(1,p1),i,drop=FALSE]
  }
  if(interaction_strength>0) {
    beta_avg = rowMeans(beta)[1:p1]
    beta_support = sort(which(abs(beta_avg) > 0))
    pairs = lapply(seq(1,round(p1/2)), function(tmp) sample(beta_support, 2, replace=FALSE))
    for(pair in pairs) {
      y = y + interaction_strength * (X[,pair[1]] * X[,pair[2]] - mean( X[,pair[1]] * X[,pair[2]]) )
    }
  }
  y = y + rnorm(n)
  return(y)
}

## Set random seed for the model
set.seed(seed.model)

## Sample the effect matrix
beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)
rowMeans(beta)

## Set random seed for the experiment
set.seed(seed)

## Sample the data
X = sample_X(n, p, p1)

#Y = sample_Y_linear(X, beta, interaction_strength=1)

set.seed(2020)
tree = generate_tree(50, p, p1)
Y = sample_Y_tree(X, beta, tree, interaction_strength=0)

# Pick a variable
j = 3
prob_treat = X[,j+p1]

res = compute.crt.ci(Y, X, j, prob_treat, seed=seed, alpha=0.1, 
                     fast=TRUE, c=50/sqrt(n), n.steps=200, roll.width=50, 
                     beta.min=-5, beta.max=5,
                     model.class="glmnet")
  
res.hist = as_tibble(res$history)
colnames(res.hist) = c("Step", "Lower", "Upper")

res.hist %>% 
  gather(Lower, Upper, key="Bound", value="Value") %>%
  ggplot(aes(x=Step, y=Value, color=Bound)) +
  geom_line() +
  geom_hline(yintercept=rowMeans(beta)[j]) +
  theme_bw()
