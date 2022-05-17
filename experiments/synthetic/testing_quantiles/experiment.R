#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_sharp.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 1000         # number of observations
p = 100          # number of variables (should be even)
signal_mean = 2  # average effect size for causal treatments
signal_std = 1   # standard deviation of effect size for causal treatments
batch <- 1
batch.size <- 1

#####################
## Input arguments ##
#####################

n = as.integer(args[1])
p = as.integer(args[2])
signal_mean = as.numeric(args[3]) / 100
signal_std = as.numeric(args[4]) / 100
batch = as.numeric(args[5])
batch_size = as.numeric(args[6])

## Parameters for P(X)
rho = 0.5
p1 = round(p/10)      # Number of treatments
p0 = p - p1           # Number of covariates

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

######################
## Experiment setup ##
######################

## Experiment header
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, rho=rho)

## Output file
out.dir = "results"
out.file = sprintf("%s/estimation_n%d_p%d_a%d-%d_b%d.txt",
                    out.dir, n, p, round(signal_mean*100), round(signal_std*100), batch)

## What should we test?
hypotheses_a = tibble(effect.quantiles = c(0.75,0.9),
                      tail = c("right", "right"))
hypotheses_b = tibble(delta = c(0,0.5*signal_mean,signal_mean,2*signal_mean))
hypotheses = crossing(hypotheses_a, hypotheses_b)


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

sample_Y = function(X, beta) {
    n = nrow(X)
    p = ncol(X)
    y = rep(0,n)
    for(i in 1:n) {
        y[i] = X[i,,drop=FALSE] %*% beta[,i,drop=FALSE]
    }
    y = y + rnorm(n)
    return(y)
}

##############
## Analysis ##
##############

analysis_test = function(X, Y, p0, j, k, tail, delta, multivariate=FALSE, reference="auto", beta=NULL) {
    ## Treatment probabilities
    prob_treat = X[,j+p1]
    ## Compute p-value
    pval = compute_pval(X, Y, j, prob_treat, k, delta, multivariate=multivariate, tail=tail, reference=reference)
    ## Check whether the hypothesis is true
    if(is.null(beta)) {
        null = NA
        tau_k = NA
        causal = NA
    } else {
        oracle = check_hypothesis(beta, j, k, delta, tail)
        null = oracle$is_null
        tau_k = oracle$tau_k
        causal = oracle$is_causal
    }
    ## Return results
    res = tibble(j=j, causal=causal, k=k, delta=delta, tail=tail, tau_k=tau_k, null=null, method="MCRT", ref=reference, multivar=multivariate, pval=pval)
    return(res)
}

experiment = function(j, seed, seed.model=2020) {

    ## Set random seed for the model
    set.seed(seed.model)

    ## Sample the effect matrix
    beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)

    ## Set random seed for the experiment
    set.seed(seed)

    ## Sample the data
    X = sample_X(n, p, p1)
    Y = sample_Y(X, beta)
    prob_treat = X[,j+p1]

    ## Define the null hypothesis
    results = tibble()
    for(h_row in 1:nrow(hypotheses)) {
        effect.quantile = hypotheses$effect.quantiles[h_row]
        k = round(effect.quantile * n)
        tail = hypotheses$tail[h_row]
        delta = hypotheses$delta[h_row]

        ## Test the null hypothesis
        for(multivariate in c(TRUE, FALSE)) {
            res = analysis_test(X, Y, p0, j, k, tail, delta, multivariate=multivariate, reference="auto", beta=beta) %>%
                mutate(seed=seed, seed.model=seed.model)
            results = rbind(results, res)
        }
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
