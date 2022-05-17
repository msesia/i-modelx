#!/usr/bin/env Rscript

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 2000         # number of observations
p = 20          # number of variables (should be even)
signal_mean = 4  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
num_inter = 2    # Interaction design
delta_int = 0.05
model.class <- "glmnet"
batch <- 1
batch.size <- 1

## Parameters for P(X)
rho = 0.5
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

## Other parameters
fdr.nominal <- 0.1
model <- "linear"

######################
## Experiment setup ##
######################

## Experiment header
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, num_inter=num_inter, delta_inter=delta_int, rho=rho)

## Output file
out.dir = "data"

#######################
## Data distribution ##
#######################


sample_interaction_design <- function(p, p1, num_inter=num_inter) {
    ## Number of non-treatment variables
    p0 = p - 2*p1
    stopifnot(p1 <= p0)
    ## Determine which interaction covariates are observed in each cluster
    return(t(sapply(1:p1, function(j) {sample(p1, num_inter)})))
}

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
    nonzero_covar = sort(sample(p0, p_causal_covar))
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

sample_Y_linear <- function(X, p1, beta_e) {
    n = nrow(X)
    p = ncol(X)
    p0 = p - 2*p1
    y = rep(0,n)
    for(i in 1:n) {
        y[i] = X[i,,drop=FALSE] %*% beta_e[,i,drop=FALSE]
    }
  y = y + rnorm(n)
  return(y)
}

calculate_effective_beta <- function(beta, Z1, idx_treat, interaction_design=interaction_design) {
    n = ncol(beta)
    beta_e = beta
    p1 = length(idx_treat)
    for(j in 1:p1) {
        idx_int = interaction_design[j,]
        for(k in idx_int) {
            beta_e[idx_treat[j],] = beta_e[idx_treat[j],] * Z1[,k]
        }
    }
    return(beta_e)
}

###########################
## Generate the data set ##
###########################

## Set random seed for the model
seed.model <- 1
set.seed(seed.model)

## Sample the effect matrix
beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)

## Sample the interaction design
interaction_design = sample_interaction_design(p, p1, num_inter=num_inter)

## Set random seed for the experiment
seed <- 1
set.seed(seed)

## Sample the data
X.full = sample_X(n, p, p1)
Z0 <- X.full[,1:p0]
Z1 <- X.full[,(p0+1):(p0+p1)]
Z <- cbind(Z0, Z1)
X <- X.full[,(p0+p1+1):p]
idx_covar <- (1):(p0+p1)
idx_treat <- (p0+p1+1):p

## Calculate effective beta
beta_e = calculate_effective_beta(beta, Z1, idx_treat, interaction_design=interaction_design)

## Sample the response
Y = sample_Y_linear(X.full, p1, beta_e)

## Save the data
df <- as_tibble(as.data.frame(cbind(Y,X,Z)))
colnames(df) <- c("Y", paste("X", 1:ncol(X), sep="_"), paste("Z", 1:ncol(Z), sep="_"))
out.file = sprintf("%s/data_n%d_p%d_a%d-%d_i%d_delta%s_%s_%s_b%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), model, model.class, batch)
df %>% write_delim(out.file, delim=" ")

## Save the ground truth
ite <- as_tibble(as.data.frame(t(beta_e)))
ite.X <- ite[,idx_treat]
ite.Z <- ite[,idx_covar]
ite <- cbind(ite.X, ite.Z)
colnames(ite) <- c(paste("X", 1:ncol(X), sep="_"), paste("Z", 1:ncol(Z), sep="_"))
out.file.ite = sprintf("%s/ite_n%d_p%d_a%d-%d_i%d_delta%s_%s_%s_b%d.txt",
                       out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), model, model.class, batch)
ite %>% write_delim(out.file.ite, delim=" ")


############################
## Generate the knockoffs ##
############################

generate.knockoff.treat <- function(X,Z) {
    p1 = ncol(X)
    U2 = matrix(runif(p1*n), n, p1)
    X.k = 0 + (U2 <= Z[,1:p1])
    return(X.k)
}
X.k <- generate.knockoff.treat(X,Z)

## Save the knockoffs
df.k <- as_tibble(as.data.frame(X.k))
colnames(df.k) <- paste("Xk", 1:ncol(X.k), sep="_")
out.file = sprintf("%s/knockoffs_n%d_p%d_a%d-%d_i%d_delta%s_%s_%s_b%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), model, model.class, batch)
df.k %>% write_delim(out.file, delim=" ")
