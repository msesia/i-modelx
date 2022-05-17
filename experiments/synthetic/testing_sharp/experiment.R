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
p = 10          # number of variables (should be even)
signal_mean = 2  # average effect size for causal treatments
signal_std = 1   # standard deviation of effect size for causal treatments
interaction_strength = 0
model <- "linear"
model.class <- "glmnet"
batch <- 1
batch.size <- 1

#####################
## Input arguments ##
#####################

n = as.integer(args[1])
p = as.integer(args[2])
signal_mean = as.numeric(args[3]) / 100
signal_std = as.numeric(args[4]) / 100
interaction_strength = as.numeric(args[5]) / 100
model <- as.character(args[6])
model.class <- as.character(args[7])
batch = as.numeric(args[8])
batch_size = as.numeric(args[9])

## Parameters for P(X)
rho = 0.8
p1 = round(p/5)       # Number of treatments
p0 = p - p1           # Number of covariates

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

######################
## Experiment setup ##
######################

## Experiment header
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, interaction_strength=interaction_strength, rho=rho)

## Output file
out.dir = "results/estimation_sharp"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_%s_%s_b%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), round(interaction_strength*100), model, model.class, batch)

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

sample_Y_tree = function(X, beta, interaction_strength=0) {
  n = nrow(X)
  p = ncol(X)
  y = rep(0,n)
  df = as_tibble(X[,seq(p1+1,p)])
  #y = predict(tree, df)
  y = rep(0,n)
  for(i in 1:n) {
    y[i] = y[i] + X[i,seq(1,p1),drop=FALSE] %*% beta[seq(1,p1),i,drop=FALSE]
  }
  #var.nonlinear <- sample(seq(p1+1,p),5,replace=FALSE)
  var.nonlinear <- c(p1+1)
  print(var.nonlinear)
  for(v in var.nonlinear) {
      y = y + 10*(X[,v]<0.05)+((X[,v]>0.5))
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

##############
## Analysis ##
##############

analysis_estimate = function(X, Y, p0, j, fast=FALSE, beta=NULL) {
    ## Treatment probabilities
    prob_treat = X[,j+p1]
    ## Compute CI
    res.crt = compute.crt.ci(Y, X, j, prob_treat, seed=seed, alpha=0.1,
                            fast=fast, c=50/sqrt(n), n.steps=200, roll.width=50, beta.min=-5, beta.max=5,
                            model.class=model.class)

    ci.crt <- res.crt$interval
    estimate.crt <- res.crt$estimate
    out.tmp <- tibble(Estimate=estimate.crt, Lower=ci.crt[1], Upper=ci.crt[2], Status=res.crt$status, Class=model.class)

    ## Check whether the hypothesis is true
    if(is.null(beta)) {
        ate = NA
    } else {
        ate = rowMeans(beta)[j]
    }
    ## Return results
    res = out.tmp %>% mutate(j=j, ate=ate, method="MCRT", fast=fast)
    return(res)
}

experiment = function(j, seed, seed.model=2020) {

    ## Set random seed for the model
    set.seed(seed.model)
   
    ## Sample the effect matrix
    beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)

    if(model=="tree") {
        set.seed(seed.model)
        tree = generate_tree(p, p, p1)
    }
    
    ## Set random seed for the experiment
    set.seed(seed)

    ## Sample the data
    X = sample_X(n, p, p1)
    if(model=="linear") {
        Y = sample_Y_linear(X, beta, interaction_strength=interaction_strength)
    } else {
        Y = sample_Y_tree(X, beta, interaction_strength=0)
    }
    prob_treat = X[,j+p1]

    ## Define the null hypothesis
    results = tibble()
    ## Test the null hypothesis
    if(model.class=="glmnet") {
        fast.list = c(TRUE,FALSE)
    } else {
        fast.list = c(FALSE)
    }
    for(fast in fast.list) {
        res = analysis_estimate(X, Y, p0, j, fast=fast, beta=beta) %>%
            mutate(seed=seed, seed.model=seed.model, Model=model)
        results = rbind(results, res)
    }

    ## Prepare output
    out = cbind(header, results)

    return(out)
}

## Run experiments
results = tibble()
for(b in 1:batch.size) {
    seed = batch.size*(batch-1) + b
    cat(sprintf("Running experiment %d of %d with seed %d...\n", b, batch.size, seed))

    for(j in 1:p1) {
        cat(sprintf("Analyzing treatment variable %d of %d......\n", j, p1))
        start_p = Sys.time()
        res = experiment(j, seed, seed.model=2020)
        end_p = Sys.time()
        print(end_p-start_p)

        ## Store results
        results = rbind(results, res)
        results %>%
            mutate_if(is.numeric, round, 4) %>%
            write_delim(out.file, delim=" ")
    }
    cat(sprintf("Written results to: %s\n", out.file))


}

cat(sprintf("All %d experiments completed.\n", batch.size))
