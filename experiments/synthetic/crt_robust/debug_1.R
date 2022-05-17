#!/usr/bin/env Rscript
options(width=200)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/methods_quantiles.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 1000         # number of observations
p = 10          # number of variables (should be even)
signal_mean = 2  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
num_inter = 1 # Interaction design
delta_int = 0.05
tail <- "left"
alpha <- 0.1
k.rel <- 10

## Set the quantile
k <- round(n*k.rel/100)

## Parameters for P(X)
rho = 0
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)
treat.threshold <- 0

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

######################
## Experiment setup ##
######################

## Output file
out.dir = "results/crtci_1"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_delta%s_%s_alpha%d_k%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), tail, round(100*alpha), k.rel)

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

## Set random seed for the model
set.seed(seed.model)

## Sample the effect matrix
beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)
idx_covar <- (p0+1):(p0+p1)
idx_treat <- (p0+p1+1):p

## Sample the interaction design
interaction_design = sample_interaction_design(p, p1, num_inter=num_inter)

## DEBUG
reference = "auto"
multivariate = FALSE
j.rel = 2

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

## Calculate effective beta
beta_e = calculate_effective_beta(beta, Z1, idx_treat, interaction_design=interaction_design)

## Sample the response
if(TRUE) {
    Y = sample_Y_linear(X.full, p1, beta_e, interaction_strength=0)
} else {
    source("../testing_sharp/nonlinear_model.R")
    set.seed(2020)
    tree = generate_tree(50, p, p1)
    Y = as.numeric(sample_Y_tree(X.full, beta_e, tree, interaction_strength=0))
}

if(FALSE) {
    ## Compute the quantile CI
    bound <- "upper"
    k <- round(0.1*n)
    check <- check_hypothesis(beta_e, j, k, 0, ifelse(bound=="lower", "right", "left"), verbose=TRUE)

    tmp = compute_interval(X.full, Y, j, prob_treat, k, multivariate=FALSE, bound=bound, reference=reference,
                           n.steps = 5000, step.size=0.5, alpha=alpha, verbose=FALSE)
}

## Compute the (regular) CI
env <- new.env()
sys.source("../utils/methods_sharp_ci.R", env)
tmp.2 <- env$compute.crt.ci(X.full, Y, j, prob_treat, k, fast=FALSE, n.steps = 1000, step.size=0.5, alpha=alpha, model.class="glmnet.inter")

tau.mean <- mean(beta_e[j,])
tau.median <- median(beta_e[j,])

tmp.2$history %>%
    gather(Lower, Upper, key="Key", value="Bound") %>%
    ggplot(aes(x=Step, y=Bound, color=Key)) +
    geom_line() +
    geom_hline(yintercept=tau.mean) +
    theme_bw()
