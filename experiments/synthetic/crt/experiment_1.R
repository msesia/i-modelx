#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_quantiles.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 500         # number of observations
p = 10          # number of variables (should be even)
signal_mean = 4  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
num_inter = 2 # Interaction design
delta_int = 0.05
tail <- "right"
delta <- 2
k.rel <- 90

#####################
## Input arguments ##
#####################
n = as.integer(args[1])
p = as.integer(args[2])
signal_mean = as.numeric(args[3]) / 100
signal_std = as.numeric(args[4]) / 100
num_inter = as.numeric(args[5])
delta_int = as.numeric(args[6]) / 100
tail <- as.character(args[7])
delta = as.numeric(args[8]) / 100
k.rel = as.numeric(args[9])

## Set the quantile
k <- round(n*k.rel/100)

## Parameters for P(X)
rho = 0.5
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

## Other parameters
##offset <- 1
##fdr.nominal <- 0.1
model <- "linear"
seed <- 2022
seed.model <- 2022
model.class <- "glmnet"
batch <- 1
batch.size <- 100


## Experiment header
tail.sign <- ifelse(tail=="right", "<=", ">=")
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, num_inter=num_inter, delta_inter=delta_int, rho=rho,
                delta=delta, k.rel=k.rel, k=k, tail.sign=tail.sign)

#######################
## Data distribution ##
#######################

sample_X <- function(n, p, p1) {
    ## Number of non-treatment variables
    p0 = p - 2*p1
    stopifnot(p1 <= p0)
    ## Generate p0 variables from a multivariate normal distribution
    mu = rep(0,p0)
    Sigma = toeplitz(rho^(0:(p0-1)))
    X1 = matrix(rnorm(n*p0),n) %*% chol(Sigma)
    X1 = pnorm(X1-1)
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
    nonzero_covar = p1 + sort(sample(p0-p1, p_causal_covar))
    nonzero = sort(c(nonzero_treat, nonzero_covar))
    # Flip signs for covariates
    signs = 2*rbinom(p,1,0.5)-1
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

## Set random seed for the model
set.seed(seed.model)

## Sample the effect matrix
beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)
idx_covar <- (p0+1):(p0+p1)
idx_treat <- (p0+p1+1):p

######################
## Experiment setup ##
######################

## Output file
out.dir = "results/crt_1"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_delta%s_%s_%d_k%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), tail, round(100*delta), k.rel)

experiment <- function(j.rel, seed, multivariate=FALSE, reference="treated") {

    ## Choose which variable to analyze
    j <- idx_treat[j.rel]

    ## Set random seed for the experiment
    set.seed(seed)

    ## Sample the data
    X.full = sample_X(n, p, p1)
    Z0 <- X.full[,1:p0]
    Z1 <- X.full[,(p0+1):(p0+p1)]
    Z <- cbind(Z0, Z1)
    X <- X.full[,(p0+p1+1):p]
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p

    ## Sample the response
    Y = sample_Y_linear(X.full, p1, beta, interaction_strength=0)

    ## Test the null hypothesis
    prob_treat <- Z0[,j.rel]
    compute_pval(X.full, Y, j, prob_treat, k, delta, multivariate=multivariate, tail=tail, reference = reference, K=1000)
}

## Run experiments
results = tibble()
for(b in 1:batch.size) {
    for(multivariate in c(FALSE, TRUE)) {
        for(reference in c("treated", "control", "auto")) {
            for(j.rel in 1:p1) {

                seed = batch.size*(batch-1) + b
                cat(sprintf("Running experiment (treatment variable %d, multivariate? %s, reference? %s) %d of %d with seed %d...\n",
                            j.rel, multivariate, reference, b, batch.size, seed))
                flush.console()

                ## Define the null hypothesis
                j <- idx_treat[j.rel]
                tmp <- check_hypothesis(beta, j, k, delta, tail, verbose=TRUE)
                hyp.str <- sprintf("Null hypothesis: %.2f = tau_(%d) %s %.3f. True? %s", tmp$tau_k, k, tail.sign, delta, tmp$is_null)

                start_p = Sys.time()

                if(FALSE) {
                    source("../utils/methods_quantiles.R")
                    experiment(j.rel, seed, multivariate=multivariate, reference="treated")
                    experiment(j.rel, seed, multivariate=multivariate, reference="control")
                }
                
                pval = experiment(j.rel, seed, multivariate=multivariate, reference=reference)
                
                end_p = Sys.time()
                print(end_p-start_p)

                cat(sprintf("Completed experiment %d of %d with seed %d.\n", b, batch.size, seed))
                res = tibble(j=j.rel, seed=seed, pval=pval, multivariate=multivariate, reference=reference, null=tmp$is_null)
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
