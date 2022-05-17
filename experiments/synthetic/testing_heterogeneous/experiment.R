#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/clustering.R")
source("../utils/mekf.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 500         # number of observations
p = 100          # number of variables (should be even)
signal_mean = 2  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
num_inter = 2 # Interaction design
delta_int = 0.05
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
num_inter = as.numeric(args[5])
delta_int = as.numeric(args[6]) / 100
model.class <- as.character(args[7])
batch = as.numeric(args[8])
batch_size = as.numeric(args[9])

## Parameters for P(X)
rho = 0.5
p1 = round(p/5)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

## Other parameters
offset <- 1
fdr.nominal <- 0.1
model <- "linear"

######################
## Experiment setup ##
######################

## Experiment header
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, num_inter=num_inter, delta_inter=delta_int, rho=rho)

## Output file
out.dir = "results/testing_heterogeneous"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_delta%s_%s_%s_b%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_inter, round(100*delta_int), model, model.class, batch)

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

sample_Y_linear <- function(X, p1, beta_e, interaction_strength=0) {
    n = nrow(X)
    p = ncol(X)
    p0 = p - 2*p1
    y = rep(0,n)
    for(i in 1:n) {
        y[i] = X[i,,drop=FALSE] %*% beta_e[,i,drop=FALSE]
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
        idx_int = interaction_design[j,]
        for(k in idx_int) {
            beta_e[idx_treat[j],] = beta_e[idx_treat[j],] * Z1[,k]
        }
    }
    return(beta_e)
}

##############
## Analysis ##
##############

knockoff.random.swap <- function(X, X.k, V.swap) {
    p1 = ncol(X)
    n = nrow(X)
    X.aug <- cbind(X, X.k)
    X.out <- X.aug
    X.out.1 <- X*(1-V.swap)+X.k*(V.swap)
    X.out.2 <- X*(V.swap)+X.k*(1-V.swap)
    X.out <- cbind(X.out.1, X.out.2)
    return(X.out)
}

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
    if(offset>1 | offset<0) {
        stop('Input offset must be between 0 or 1')
    }
    ts = sort(c(0, abs(W)))
    ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
    ok = which(ratio <= fdr)
    ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

define_hypotheses <- function(vars.split, beta_e, Z, X) {
    p.x <- ncol(X)
    p0 <- ncol(Z) - p.x
    n <- nrow(X)
    ## List of individual clusters (environments) for each variable
    env.list <- lapply(1:p1, function(j) { partition_covariates(Z, vars.split[[j]])})
    V.list <- lapply(1:p.x, function(j) {
        num_env = max(env.list[[j]])
        w.j <- rep(j, num_env)
        return(w.j)
    })
    E.list <- lapply(1:p.x, function(j) {
        num_env = max(env.list[[j]])
        w.j <- 1:num_env
        return(w.j)
    })

    hyp <- tibble(Treatment=unlist(V.list), Partition=unlist(E.list))
    hyp$Splits <- sapply(1:nrow(hyp), function(i) {length(vars.split[[hyp$Treatment[i]]])})
    hyp$Label <- sapply(1:nrow(hyp), function(i) {
        if(length(vars.split[[hyp$Treatment[i]]])>0) {
            tmp.1 <- paste("T", vars.split[[hyp$Treatment[i]]], sep="")
            tmp.2 <- binaryLogic::as.binary(hyp$Partition[i]-1,n=hyp$Splits[i])
            tmp.3 <- paste(tmp.1, tmp.2,sep=":")
            out <- paste(tmp.3, collapse=",")
        } else {
            out <- "None"
        }
        return(out)
    })
    hyp <- hyp %>% select(-Splits)
    ## Check whether hypotheses are null
    hyp$Null <- NA
    hyp$Causal.prop <- NA
    hyp$Beta.sd <- NA
    for(i in 1:nrow(hyp)) {
        j <- hyp$Treatment[i]
        idx <- which(env.list[[j]]==hyp$Partition[i])
        hyp$Null[i] = all(beta_e[p0+p1+j,idx]==0)
        hyp$Causal.prop[i] = mean(beta_e[p0+p1+j,idx]!=0)
        hyp$Beta.sd[i] = sd(beta_e[p0+p1+j,idx])
    }

    return(hyp)
}

evaluate_results <- function(discoveries, hypotheses) {
    ## Evaluate FDP
    if(nrow(discoveries)>0) {
        discoveries.1 <- left_join(discoveries, hypotheses, by = c("Treatment", "Partition", "Label"))
        res.tmp <- discoveries.1 %>% summarise(FDP=mean(Null), Causal.prop=mean(Causal.prop), Beta.sd=mean(Beta.sd))
    } else {
        res.tmp <- tibble(FDP=0, Causal.prop=NA, Beta.sd=NA)
    }
    ## Evaluate power
    if(nrow(discoveries)>0) {
        discoveries.2 <- left_join(hypotheses, mutate(discoveries, Discovered=TRUE), by = c("Treatment", "Partition", "Label")) %>%
            mutate(Discovered=ifelse(is.na(Discovered), FALSE, Discovered))
        res.tmp$Power <- discoveries.2 %>% filter(!Null) %>% summarise(Power=mean(Discovered)) %>% as.numeric()
        res.tmp$True <- discoveries.2 %>% filter(!Null) %>% summarise(Power=sum(Discovered)) %>% as.numeric()
        res.tmp$Discoveries <- discoveries.2 %>% summarise(Power=sum(Discovered)) %>% as.numeric()
    } else {
        res.tmp$Power <- 0
        res.tmp$True <- 0
        res.tmp$Discoveries <- 0
    }
    return(res.tmp)
}

analysis_knockoffs <- function(Z, X, Y, q=0.1, offset=0, naive=FALSE, split=FALSE, vanilla=FALSE, method.cluster="Lasso") {
    n = nrow(X)
    p1 = ncol(X)
    Z0 <- Z[,1:p0]
    Z1 <- Z[,(p0+1):(p0+p1)]

    cross.prior=TRUE
    if(naive) cross.prior = FALSE
    if(split) cross.prior = FALSE
    if(vanilla) cross.prior = FALSE

    random.swap=TRUE
    if(naive) random.swap = FALSE
    if(vanilla) random.swap = FALSE

    
    generate.knockoff.treat <- function(X) {
        p1 = ncol(X)
        U2 = matrix(runif(p1*n), n, p1)
        X.k = 0 + (U2 <= Z[,1:p1])
        return(X.k)
    }
    X.k <- generate.knockoff.treat(X)

    ## Randomly swap variables and knockoffs
    if(split) {
        fold.1 <- sort(sample(n,round(n/2),replace=FALSE))
        fold.2 <- setdiff(1:n, fold.1)
        V.swap.1 <- matrix(rep(0,n*p1), n)
        V.swap.2 <- matrix(rbinom(n*p1,1,1/2), n)
    } else {
        if(naive) {
            V.swap <- matrix(rep(0,n*p1), n)
        } else {
            V.swap <- matrix(rbinom(n*p1,1,1/2), n)
        }
        fold.1 <- 1:n
        fold.2 <- 1:n
        V.swap.1 <- V.swap
        V.swap.2 <- V.swap
    }
    X.swap = knockoff.random.swap(X, X.k, V.swap.1)

    ## Clustering HTE using knockoff-augmented data
    if(vanilla) {
        env.list <- lapply(1:p1, function(j) { rep(1,n) })
        vars.split <- lapply(1:p1, function(j) { c() })
    } else {
        cat(sprintf("Searching for interactions...\n"))
        flush.console()
        if(method.cluster=="BART") {
            cluster.res <- cluster_covariates_BART(Z0[fold.1,], Z1[fold.1,], X.swap[fold.1,], Y[fold.1], delta_int=delta_int)
        } else {
            cluster.res <- cluster_covariates_lasso(Z0[fold.1,], Z1[fold.1,], X.swap[fold.1,], Y[fold.1])
        }
        ## List of binary variable splits for each treatment
        vars.split <- lapply(1:p1, function(j) { filter(cluster.res,treatment==j)$covariate})
        ## List of individual clusters (environments) for each variable
        env.list <- lapply(1:p1, function(j) { partition_covariates(Z[fold.2,], vars.split[[j]])})
        cat(sprintf("Done searching for interactions.\n"))
        flush.console()
    }
    
    ## Compute multi-environment knockoff statistics
    cat(sprintf("Computing test statistics...\n"))
    flush.console()
    stats <- compute_stats_by_group(Z[fold.2,], X[fold.2,], X.k[fold.2,], Y[fold.2], V.swap.2[fold.2,], env.list, family = "gaussian",
                                    cross.prior=cross.prior, random.swap=random.swap, verbose=TRUE, dfmax=500)
    cat(sprintf("Done computing test statistics.\n"))
    flush.console()
    
    ## Vectorize the statistics
    W.vec <- unlist(stats$W)
    V.vec <- unlist(stats$V)
    E.vec <- unlist(stats$E)

    ## Apply the knockoff filter
    thres = knockoff::knockoff.threshold(W.vec, fdr=fdr.nominal, offset=offset)
    S.vec <- which(W.vec >= thres)

    ## Make the results intellegible
    if(length(S.vec)>0) {
        discoveries <- tibble(Treatment=V.vec[S.vec], Partition=E.vec[S.vec])
        discoveries$Splits <- sapply(1:nrow(discoveries), function(i) {length(vars.split[[discoveries$Treatment[i]]])})
        discoveries$Label <- sapply(1:nrow(discoveries), function(i) {
            if(length(vars.split[[discoveries$Treatment[i]]])>0) {
                tmp.1 <- paste("T", vars.split[[discoveries$Treatment[i]]], sep="")
                tmp.2 <- binaryLogic::as.binary(discoveries$Partition[i]-1,n=discoveries$Splits[i])
                tmp.3 <- paste(tmp.1, tmp.2,sep=":")
                out <- paste(tmp.3, collapse=",")
            } else {
                out <- "None"
            }
            return(out)
        })
        discoveries <- discoveries %>% select(-Splits)
    } else {
        discoveries <- tibble(Treatment=NA, Partition=NA, Label=NA) %>% head(0)
    }

    cat(sprintf("List of %d discoveries:\n", nrow(discoveries)))
    flush.console()
    print(discoveries)
    
    ## Output results
    out <- c()
    out$discoveries <- discoveries
    out$splits <- vars.split
    return(out)
}

experiment <- function(seed, seed.model=2022) {

    ## Set random seed for the model
    set.seed(seed.model)

    ## Sample the effect matrix
    beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)

    ## Sample the interaction design
    interaction_design = sample_interaction_design(p, p1, num_inter=num_inter)

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

    ## Calculate effective beta
    beta_e = calculate_effective_beta(beta, Z1, idx_treat, interaction_design=interaction_design)

    ## Sample the response
    Y = sample_Y_linear(X.full, p1, beta_e, interaction_strength=0)

    #############################
    ## EXACT knockoff analysis ##
    #############################
    cat(sprintf("Applying EXACT method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.exact <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset)
    discoveries.exact <- out.exact$discoveries
    ## Define hypotheses
    vars.split <- out.exact$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    ## Assess results
    res.exact <- evaluate_results(discoveries.exact, hypotheses) %>% mutate(Method="Exact", q=fdr.nominal, offset=offset)
    cat(sprintf("Done with EXACT method.\n"))
    flush.console()

    #############################
    ## NAIVE knockoff analysis ##
    #############################
    cat(sprintf("Applying NAIVE method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.naive <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, naive=TRUE)
    discoveries.naive <- out.naive$discoveries
    ## Define hypotheses
    vars.split <- out.naive$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    ## Assess results
    res.naive <- evaluate_results(discoveries.naive, hypotheses) %>% mutate(Method="Naive", q=fdr.nominal, offset=offset)
    cat(sprintf("Done with NAIVE method.\n"))
    flush.console()

    #############################
    ## SPLIT knockoff analysis ##
    #############################
    cat(sprintf("Applying SPLIT method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.split <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, split=TRUE)
    discoveries.split <- out.split$discoveries
    ## Define hypotheses
    vars.split <- out.split$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    ## Assess results
    res.split <- evaluate_results(discoveries.split, hypotheses) %>% mutate(Method="Split", q=fdr.nominal, offset=offset)
    cat(sprintf("Done with SPLIT method.\n"))
    flush.console()

    ###############################
    ## VANILLA knockoff analysis ##
    ###############################
    cat(sprintf("Applying VANILLA method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.vanilla <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, vanilla=TRUE)
    discoveries.vanilla <- out.vanilla$discoveries
    ## Define hypotheses
    vars.split <- out.vanilla$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    ## Assess results
    res.vanilla <- evaluate_results(discoveries.vanilla, hypotheses) %>% mutate(Method="Vanilla", q=fdr.nominal, offset=offset)
    cat(sprintf("Done with VANILLA method.\n"))
    flush.console()

    ## Combine results
    out <- rbind(res.exact, res.naive, res.split, res.vanilla)
    ##out <- rbind(res.split)

    return(out)
}

## Run experiments
results = tibble()
for(b in 1:batch.size) {
    seed = batch.size*(batch-1) + b
    cat(sprintf("Running experiment %d of %d with seed %d...\n", b, batch.size, seed))
    flush.console()
    
    start_p = Sys.time()
    res = experiment(seed, seed.model=2022)
    end_p = Sys.time()
    print(end_p-start_p)

    cat(sprintf("Completed experiment %d of %d with seed %d.\n", b, batch.size, seed))
    print(res)
    flush.console()
    
    ## Store results
    results = rbind(results, res)
    cbind(results,header) %>%
        mutate_if(is.numeric, round, 4) %>%
        write_delim(out.file, delim=" ")

    cat(sprintf("Written results to: %s\n", out.file))
    flush.console()
}

cat(sprintf("All %d experiments completed.\n", batch.size))
flush.console()
