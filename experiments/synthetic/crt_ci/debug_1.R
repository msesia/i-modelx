#!/usr/bin/env Rscript

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_quantiles.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 2000         # number of observations
p = 20          # number of variables (should be even)
signal_mean = 4  # average effect size for causal treatments
signal_std = 1   # standard deviation of effect size for causal treatments
num_inter = 2 # Interaction design
delta_int = 0.05
tail <- "left"
alpha <- 0.1
k.rel <- 10

## Set the quantile
k <- round(n*k.rel/100)

## Parameters for P(X)
rho = 0.5
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)
treat.threshold <- 1

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
                k.rel=k.rel, k=k, tail.sign=tail.sign, bound=bound, alpha=alpha, treat.threshold=treat.threshold, prop.treat=prop.treat)

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
out.dir = "results/crtci_1"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_delta%s_%s_alpha%d_k%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), tail, round(100*alpha), k.rel)

if(FALSE) {
    j.rel <- 3
    b <- 1
    seed <- 1
    reference = "auto"
    multivariate = FALSE
}

experiment <- function(j.rel, seed, alpha=0.1, multivariate=FALSE, reference="treated", bound="lower") {


    if(FALSE) {
        source("../utils/methods_quantiles.R")
        k <- round(0.1*n)
        check <- check_hypothesis(beta, j, k, 0, "right", verbose=TRUE)
        tmp = compute_interval(X.full, Y, j, prob_treat, k, multivariate=TRUE, bound="lower", reference="auto",
                               n.steps = 1000, step.size=1, alpha=alpha, verbose=FALSE)
        plot(tmp$history$step, tmp$history$bound); abline(h=check$tau_k)

        k <- round(0.1*n)        
        tmp = compute_interval(X.full, Y, j, prob_treat, k, multivariate=TRUE, bound="upper", reference="auto",
                               n.steps = 1000, step.size=1, alpha=alpha, verbose=FALSE)
        plot(tmp$history$step, tmp$history$bound)

    }

    if(FALSE) {
        source("../utils/methods_quantiles.R")
        tail <- "left"
        k <- round(0.9*n)
        delta <- 20
        tmp <- check_hypothesis(beta, j, k, delta, tail, verbose=TRUE)
        compute_pval(X.full, Y, j, prob_treat, k, delta, multivariate=FALSE, tail=tail, reference="control", K=1000)
        compute_pval(X.full, Y, j, prob_treat, k, delta, multivariate=FALSE, tail=tail, reference="treatment", K=1000)

        source("../utils/methods_quantiles.R")
        tail <- "right"
        k <- round(0.1*n)
        delta <- 10
        tmp <- check_hypothesis(beta, j, k, delta, tail, verbose=TRUE)
        compute_pval(X.full, Y, j, prob_treat, k, delta, multivariate=FALSE, tail=tail, reference="control", K=1000)
        compute_pval(X.full, Y, j, prob_treat, k, delta, multivariate=FALSE, tail=tail, reference="treatment", K=1000)
}
    
    bound.value = tmp$bound
    return(bound.value)
}

reference = "auto"
multivariate = FALSE
j.rel = 3

## Choose which variable to analyze
j <- idx_treat[j.rel]

## Set random seed for the experiment
#set.seed(seed)

## Sample the data
X.full = sample_X(n, p, p1, treat.threshold=treat.threshold)
Z0 <- X.full[,1:p0]
Z1 <- X.full[,(p0+1):(p0+p1)]
Z <- cbind(Z0, Z1)
X <- X.full[,(p0+p1+1):p]
idx_covar <- (p0+1):(p0+p1)
idx_treat <- (p0+p1+1):p
prob_treat <- Z0[,j.rel]

## Sample the response
Y = sample_Y_linear(X.full, p1, beta, interaction_strength=0)

## Compute the CI 
bound <- "upper"
k <- round(0.1*n)
check <- check_hypothesis(beta, j, k, 0, ifelse(bound=="lower", "right", "left"), verbose=TRUE)

tmp = compute_interval(X.full, Y, j, prob_treat, k, multivariate=FALSE, bound=bound, reference=reference,
                       n.steps = 5000, step.size=0.5, alpha=alpha, verbose=FALSE)


plot(tmp$history$step, tmp$history$bound, ylim=c(-20,40)); abline(h=check$tau_k)
