#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(width=200, error=function()traceback(2))

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_quantiles.R")
env <- new.env()
sys.source("../utils/methods_sharp_ci.R", env)

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 10         # number of observations
signal_mean = 2  # average effect size for causal treatments
signal_std = 1   # standard deviation of effect size for causal treatments
num_inter = 2 # Interaction design
model.class <- "glmnet.fast"
ci.type <- "sharp"
tail <- "NA"
k.rel <- 50
batch <- 1

#####################
## Input arguments ##
#####################
if(TRUE) {
    n = as.integer(args[1])
    signal_mean = as.numeric(args[2]) / 100
    signal_std = as.numeric(args[3]) / 100
    num_inter = as.numeric(args[4])
    model.class = as.character(args[5])
    ci.type = as.character(args[6])
    tail <- as.character(args[7])
    k.rel = as.numeric(args[8])
    batch = as.integer(args[9])
}

## Set the quantile
k <- round(n*k.rel/100)

## Parameters for P(X)
p = 10                  # number of variables (should be even)
rho = 0
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)
treat.threshold = 1

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

## Other parameters
alpha <- 0.1
model <- "linear"
seed <- 2022
seed.model <- 2022
delta_int = 0.05
batch.size <- 10

#######################
## Data distribution ##
#######################

estimate_prop_treat <- function(treat.threshold, n=1000000) {
    mu = rep(0,p0)
    Sigma = toeplitz(rho^(0:(p0-1)))
    X1 = matrix(rnorm(n*p0),n) %*% chol(Sigma)
    X1 = pnorm(X1-treat.threshold)
    round(mean(colMeans(X1)),2)
}
prop.treat = estimate_prop_treat(treat.threshold)

## Experiment header
tail.sign <- ifelse(tail=="right", "<=", ">=")
bound = ifelse(tail=="right", "lower", "upper")
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, num_inter=num_inter, delta_inter=delta_int, rho=rho,
                k.rel=k.rel, k=k, tail.sign=tail.sign, alpha=alpha, treat.threshold=treat.threshold, prop.treat=prop.treat,
                model.class=model.class, ci.type=ci.type, batch=batch)

sample_interaction_design <- function(p, p1, num_inter=num_inter) {
    if(num_inter==0) return(matrix(0,0,0))
    ## Number of non-treatment variables
    p0 = p - 2*p1
    stopifnot(p1 <= p0)
    ## Determine which interaction covariates are observed in each cluster
    out <- t(sapply(1:p1, function(j) {sample(p1, num_inter)}))
    out <- matrix(out, p1, num_inter)
    return(out)
}

sample_X <- function(n, p, p1, treat.threshold=1) {
    ## Number of non-treatment variables
    p0 = p - 2*p1
    stopifnot(p1 <= p0)
    ## Generate p0 variables from a multivariate normal distribution
    mu = rep(0,p0)
    Sigma = toeplitz(rho^(0:(p0-1)))
    X1 = matrix(rnorm(n*p0),n) %*% chol(Sigma)
    X1 = pnorm(X1-treat.threshold)
    ## Generate p1 variables from Bernoulli conditional on the first p1 of the p0 variables
    U2 = matrix(runif(p1*n), n, p1)
    X2 = (U2 <= X1[,1:p1])
    ## Generate interaction covariates
    X3 = matrix(rep(0, n*p1), n, p1)
    for(j in 1:p1) {
        X3[,j] = rbinom(n, 1, 0.5)
    }
    #cluster_means = matrix(rbeta(K*p1, 0.5, 0.5), K, p1)
    # Combine the variables (p0 covariates, p1 interaction covariates, p1 treatments)
    X = cbind(X1, X3, X2)
    return(X)
}

sample_effects <- function(p, p1, p_causal, p_causal_covar, signal_mean, signal_std) {
    p0 = p - 2*p1
    nonzero_treat = p0+p1+sort(sample(p1, p_causal))
    nonzero_covar = c() #p1 + sort(sample(p0-p1, p_causal_covar))
    nonzero = sort(c(nonzero_treat, nonzero_covar))
    # Flip signs for covariates
    signs = 1 #2*rbinom(p,1,0.5)-1
    beta_center = signal_mean * (1:p %in% nonzero) * signs
    beta_center = matrix(rep(beta_center,each=n), n, p)
    beta_matrix = beta_center + matrix(rnorm(n*p, mean=0, sd=signal_std), n, p)
    beta_matrix[,-nonzero] = 0
    beta_matrix = t(beta_matrix)
    beta_matrix[(p0+1):(p0+p1),] = 0
    return(beta_matrix)
}

sample_Y_linear <- function(X, p1, beta, interaction_strength=0) {
    n = nrow(X)
    p = ncol(X)
    p0 = p - 2*p1
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

sample_Y_tree <- function(X, beta, tree, interaction_strength=0) {
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

calculate_effective_beta <- function(beta, Z1, idx_treat, interaction_design=interaction_design) {
    n = ncol(beta)
    beta_e = beta
    p1 = length(idx_treat)
    for(j in 1:p1) {
        if(nrow(interaction_design)>0) {
            idx_int = interaction_design[j,]
            for(k in idx_int) {
                beta_e[idx_treat[j],] = beta_e[idx_treat[j],] * Z1[,k]
            }
        }
    }
    return(beta_e)
}

######################
## Experiment setup ##
######################

## Set random seed for the model
set.seed(seed.model)

## Sample the effect matrix
beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)
idx_covar <- (p0+1):(p0+p1)
idx_treat <- (p0+p1+1):p

## Sample the interaction design
interaction_design = sample_interaction_design(p, p1, num_inter=num_inter)

## Output file
out.dir = "results/crt_robust_1"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_%s_%s_t%s_k%d_batch%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, model.class, ci.type, tail, k.rel, batch)

experiment <- function(j.rel, seed, alpha=0.1, multivariate=FALSE, reference="treated", bound="lower", ci.type="sharp") {
    ## Set random seed for the experiment
    set.seed(seed)

    ## Sample the data
    X.full = sample_X(n, p, p1, treat.threshold=treat.threshold)
    Z0 <- X.full[,1:p0]
    Z1 <- X.full[,(p0+1):(p0+p1)]
    Z <- cbind(Z0, Z1)
    X <- X.full[,(p0+p1+1):p]
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p
    prob_treat <- Z0[,j.rel]

    ## Calculate effective beta
    beta_e = calculate_effective_beta(beta, Z1, idx_treat, interaction_design=interaction_design)

    ## Sample the response
    Y = sample_Y_linear(X.full, p1, beta_e, interaction_strength=0)

    ## Choose which variable to analyze
    j <- idx_treat[j.rel]

    ## Compute the sharp CI
    if(ci.type=="sharp") {
        if(model.class=="glmnet.fast") {
            tmp <- env$compute.crt.ci(X.full, Y, j, prob_treat, k, fast=TRUE, n.steps = 1000, step.size=0.5, alpha=alpha, model.class="glmnet")
        } else {
            tmp <- env$compute.crt.ci(X.full, Y, j, prob_treat, k, fast=FALSE, n.steps = 1000, step.size=0.5, alpha=alpha, model.class=model.class)
        }
        tmp.check <- check_hypothesis(beta_e, j, round(0.5*n), 0, "lower", verbose=FALSE)
        out <- c()
        out$bound <- c("Lower", "Upper")
        out$bound.value <- tmp$interval
        out$tau_k <- tmp.check$tau_k
        out$is_null <- all(beta_e[j,]==0)
    } else {
        tmp.check <- check_hypothesis(beta_e, j, k, 0, tail, verbose=FALSE)
        tmp = compute_interval(X.full, Y, j, prob_treat, k, multivariate=multivariate, bound=bound,
                               reference=reference, n.steps = 2500, step.size=0.5, alpha=alpha, verbose=FALSE)        
        out <- c()
        out$bound <- c(bound)
        out$bound.value <- tmp$bound
        out$tau_k <- tmp.check$tau_k
        out$is_null <- tmp.check$is_null

    }
    return(out)
}


## Run experiments
results = tibble()
for(b in 1:batch.size) {
    for(multivariate in c(FALSE, TRUE)) {
        for(reference in c("auto")) {
            for(j.rel in 1:p1) {

                seed = batch.size*(batch-1) + b
                cat(sprintf("Running experiment (treatment variable %d, multivariate? %s, reference? %s) %d of %d with seed %d...\n",
                            j.rel, multivariate, reference, b, batch.size, seed))
                flush.console()

                ## Define the variable to analyze
                j <- idx_treat[j.rel]

                start_p = Sys.time()

                ci.res = experiment(j.rel, seed, alpha=alpha, multivariate=multivariate, reference=reference, bound=bound, ci.type=ci.type)

                end_p = Sys.time()
                print(end_p-start_p)

                cat(sprintf("Completed experiment %d of %d with seed %d.\n", b, batch.size, seed))
                res = tibble(j=j.rel, seed=seed,
                             bound=ci.res$bound, bound.value=ci.res$bound.value,
                             tau_k=ci.res$tau_k, multivariate=multivariate, reference=reference, null=ci.res$is_null)
                flush.console()

                ## Store results
                results = rbind(results, res)
                cbind(results,header) %>%
                    mutate_if(is.numeric, round, 4) %>%
                    write_delim(out.file, delim=" ")

                cat(sprintf("Written results to: %s\n", out.file))
                flush.console()
            }
        }
    }
}

cat(sprintf("All %d experiments completed.\n", batch.size))
flush.console()
