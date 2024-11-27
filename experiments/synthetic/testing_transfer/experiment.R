#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####################
## Load libraries ##
####################

suppressMessages(library(tidyverse))
source("../utils/clustering_ite.R")
source("../utils/mekf.R")
source("../utils/mekf_original.R")

###########################
## Simulation parameters ##
###########################

## Problem parameters
n = 1000         # number of observations
p = 160          # number of variables (should be even)
signal_mean = 4  # average effect size for causal treatments
signal_std = 0   # standard deviation of effect size for causal treatments
num_cluster = 4 # Interaction design
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
num_cluster = as.numeric(args[5])
delta_int = as.numeric(args[6]) / 100
model.class <- as.character(args[7])
batch = as.numeric(args[8])
batch_size = as.numeric(args[9])

## Parameters for P(X)
rho = 0.5
p1 = round(p/4)         # Number of treatments
p0 = p - 2*p1           # Number of covariates (constraint: 2*p1 < p0)

## Parameters for P(Y|X)
p_causal = round(p1/2)        # Number of causal treatments
p_causal_covar = round(p0/2)  # Number of causal covariates

## Other parameters
offset <- 1
fdr.nominal <- 0.1
model <- "linear"
cluster.method <- "Lasso"

## DEBUG
if(FALSE) {
    b <- 1
    seed = batch.size*(batch-1) + b
    seed.model = 2022
    naive = FALSE
    split = FALSE
    KFglobal = FALSE
    method.cluster = "BART"
}

######################
## Experiment setup ##
######################

## Experiment header
header = tibble(n=n, p=p, p_causal=p_causal, p_causal_covar=p_causal_covar,
                signal_mean=signal_mean, signal_std=signal_std, num_cluster=num_cluster, delta_inter=delta_int, rho=rho)

## Output file
out.dir = "results/testing_transfer"
out.file = sprintf("%s/n%d_p%d_a%d-%d_i%d_delta%s_%s_%s_b%d.txt",
                   out.dir, n, p, round(signal_mean*100), round(signal_std*100), num_cluster, round(100*delta_int), model, model.class, batch)

#######################
## Data distribution ##
#######################


sample_interaction_design <- function(p, p1, num_inter=2) {
    ## Number of non-treatment variables
    p0 = p - 2*p1
    stopifnot(p1 <= p0)
    ## Determine which interaction covariates are observed in each cluster
    out <- lapply(1:p1, function(j) {
        if(j%%2==0) {
            sample(p1, num_inter)
        } else {
            c()
        }
    })
    return(out)
}

sample_X <- function(n, p, p1, test=FALSE) {
    ## This will produce:
    ## p0 covariates (the first p1 of which p1 determine the propensity scores)
    ## p1 interaction covariates
    ## p1 treatments
    ## Note: p = p0+p1
    ##
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
    ## Generate interaction covariates (zeros at test time)
    X3 = matrix(rep(0, n*p1), n, p1)
    if(test==FALSE) {
        for(j in 1:p1) {
            X3[,j] = rbinom(n, 1, 0.5)
        }
    }
    #cluster_means = matrix(rbeta(K*p1, 0.5, 0.5), K, p1)
    # Combine the variables (p0 covariates, p1 interaction covariates, p1 treatments)
    X = cbind(X1, X3, X2)
    return(X)
}

sample_effects <- function(p, p1, p_causal, p_causal_covar, signal_mean, signal_std) {
    ## Structure of X:
    ## p0 covariates (the first p1 of which p1 determine the propensity scores)
    ## p1 interaction covariates (no effects)
    ## p1 treatments (of which some are causal)
    ## Note: p = p0+p1
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
        idx_int = interaction_design[[j]]
        if(!is.null(idx_int)) {
            for(k in idx_int) {
                beta_e[idx_treat[j],] = beta_e[idx_treat[j],] * Z1[,k]
            }
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

define_hypotheses <- function(env.list, beta_e, Z, X) {
    p.x <- ncol(X)
    p0 <- ncol(Z) - p.x
    n <- nrow(X)
    ## List of individual clusters (environments) for each variable
    V.list <- 1:p.x

    hyp <- tibble(Treatment=unlist(V.list))
    ## Check whether hypotheses are null
    hyp$Null <- NA
    hyp$Causal.prop <- NA
    hyp$Beta.sd <- NA
    for(i in 1:nrow(hyp)) {
        j <- hyp$Treatment[i]
        idx <- 1:n
        hyp$Null[i] = all(beta_e[p0+p1+j,idx]==0)
        hyp$Causal.prop[i] = mean(beta_e[p0+p1+j,idx]!=0)
        hyp$Beta.sd[i] = sd(beta_e[p0+p1+j,idx])
    }

    return(hyp)
}

evaluate_results_grouped <- function(discoveries, hypotheses) {
    if(nrow(discoveries)>0) {
        r.max <- num_cluster
        significance.methods <- unique(discoveries$significance)
        tmp <- lapply(significance.methods, function(s.method) {
            tmp <- lapply(1:r.max, function(r.idx) {
                evaluate_results(filter(discoveries,r==r.idx,significance==s.method), hypotheses) %>%
                    mutate(r=r.idx, significance=s.method)
            })
            tmp <- do.call("rbind", tmp)
        })
        tmp <- do.call("rbind", tmp)
    } else {
        tmp <- tibble()
    }
    return(tmp)
}

evaluate_results <- function(discoveries, hypotheses) {
    ## Evaluate FDP
    if(nrow(discoveries)>0) {
        discoveries.1 <- left_join(discoveries, hypotheses, by = c("Treatment"))
        res.tmp <- discoveries.1 %>% summarise(FDP=mean(Null), Causal.prop=mean(Causal.prop), Beta.sd=mean(Beta.sd)) %>%
            mutate(FDP=ifelse(is.na(FDP),0,FDP))
    } else {
        res.tmp <- tibble(FDP=0, Causal.prop=NA, Beta.sd=NA)
    }
    ## Evaluate power
    if(nrow(discoveries)>0) {
        discoveries.2 <- left_join(hypotheses, mutate(discoveries, Discovered=TRUE), by = c("Treatment")) %>%
            mutate(Discovered=ifelse(is.na(Discovered), FALSE, Discovered))
        res.tmp$Power <- discoveries.2 %>% filter(!Null) %>% summarise(Power=mean(Discovered)) %>% as.numeric()
        res.tmp$Power <- ifelse(is.nan(res.tmp$Power), 0, res.tmp$Power)
        res.tmp$True <- discoveries.2 %>% filter(!Null) %>% summarise(Power=sum(Discovered)) %>% as.numeric()
        res.tmp$Discoveries <- discoveries.2 %>% summarise(Power=sum(Discovered)) %>% as.numeric()
    } else {
        res.tmp$Power <- 0
        res.tmp$True <- 0
        res.tmp$Discoveries <- 0
    }
    return(res.tmp)
}

analysis_knockoffs <- function(Z, X, Y, q=0.1, offset=0, cluster.method="Lasso", naive=FALSE, split=FALSE, KFglobal=FALSE) {
    n = nrow(X)
    p1 = ncol(X)
    Z0 <- Z[,1:p0]
    Z1 <- Z[,(p0+1):(p0+p1)]

    cross.prior=TRUE
    if(naive) cross.prior = FALSE
    if(split) cross.prior = FALSE
    if(KFglobal) cross.prior = FALSE

    random.swap=TRUE
    if(naive) random.swap = FALSE
    if(KFglobal) random.swap = FALSE

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
    if(KFglobal) {
        env.list <- lapply(1:p1, function(j) { rep(1,n) })
    } else {
        cat(sprintf("Clustering covariate space...\n"))
        flush.console()
        if(cluster.method=="Lasso") {
            env.list <- cluster_covariates_lasso(Z0, Z1, X.swap, Y, fold.1, fold.2, num_cluster=num_cluster)
        } else {
            env.list <- lapply(1:p1, function(j) {
                cat(sprintf("Clustering covariate space for treatment %d...\n", j))
                flush.console()
                cluster_covariates_ite(j, Z0, Z1, X.swap, Y, fold.1, fold.2, num_cluster=num_cluster)
            })
        }
        cat(sprintf("Done clustering covariate space.\n"))
        flush.console()
    }

    ## Compute multi-environment knockoff statistics
    cat(sprintf("Computing test statistics...\n"))
    flush.console()
    stats <- compute_stats_by_group(Z[fold.2,], X[fold.2,], X.k[fold.2,], Y[fold.2], V.swap.2[fold.2,], env.list, family = "gaussian",
                                    cross.prior=cross.prior, random.swap=random.swap, verbose=TRUE, dfmax=500, nfolds=5)
    cat(sprintf("Done computing test statistics.\n"))
    flush.console()

    ## Prepare stats for multi-environment knockoff filter
    if(KFglobal) {
        W.matrix = unlist(stats$W)
    } else {
        W.stats <- lapply(1:length(stats$W), function(j) {
            w <- stats$W[[j]]
            if(length(w) < num_cluster) {
                w <- rep(w, each=num_cluster/length(w))
            }
            return(w)
        })
        W.matrix = t(do.call("cbind", W.stats))
    }
    
    ## Apply the multi-environment knockoff filter
    discoveries <- tibble()
    if(KFglobal) {
        W.vec <- as.numeric(W.matrix)
        thres <- knockoff::knockoff.threshold(W.vec, fdr=fdr.nominal, offset=offset)
        S <- which(W.vec >= thres)
        if(length(S)==0) S <- c(-1)
        discoveries.tmp <- tibble(Treatment=S, r=1, q=fdr.nominal, significance="KF")
        discoveries <- rbind(discoveries, discoveries.tmp)
    } else {
        for(r.pc in 1:num_cluster) {
            S.1 <- partial_conjunction(W.matrix, r.pc, q = fdr.nominal, method = "accumulation", randomness = 3)            
            if(length(S.1)==0) S.1 <- c(-1)
            discoveries.tmp <- tibble(Treatment=S.1, r=r.pc, q=fdr.nominal, significance="accumulation")
            discoveries <- rbind(discoveries, discoveries.tmp)
            S.2 <- partial_conjunction(W.matrix, r.pc, q = fdr.nominal, method = "seqstep", c=0.4, randomness = 3, offset=offset)
            if(length(S.2)==0) S.2 <- c(-1)
            discoveries.tmp <- tibble(Treatment=S.2, r=r.pc, q=fdr.nominal, significance="seqstep-0.4")
            discoveries <- rbind(discoveries, discoveries.tmp)
            S.3 <- partial_conjunction(W.matrix, r.pc, q = fdr.nominal, method = "seqstep", c=0.5, randomness = 3, offset=offset)
            if(length(S.3)==0) S.3 <- c(-1)
            discoveries.tmp <- tibble(Treatment=S.3, r=r.pc, q=fdr.nominal, significance="seqstep-0.5")
            discoveries <- rbind(discoveries, discoveries.tmp)
            S.4 <- partial_conjunction(W.matrix, r.pc, q = fdr.nominal, method = "seqstep", c=0.6, randomness = 3, offset=offset) 
            if(length(S.4)==0) S.4 <- c(-1)
            discoveries.tmp <- tibble(Treatment=S.4, r=r.pc, q=fdr.nominal, significance="seqstep-0.6")
            discoveries <- rbind(discoveries, discoveries.tmp)
        }
    }

    cat(sprintf("List of %d discoveries:\n", nrow(discoveries)))
    flush.console()
    print(discoveries)

    ## Output results
    out <- c()
    out$discoveries <- discoveries
    out$clusters <- env.list
    return(out)
}

experiment <- function(seed, seed.model=2022) {

    ## Set random seed for the model
    set.seed(seed.model)

    ## Sample the effect matrix
    beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)

    ## Sample the interaction design
    interaction_design = sample_interaction_design(p, p1, num_inter=2)

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

    ## Sample the test data (with covariate shift)
    X.full.test = sample_X(n, p, p1, test=TRUE)
    Z0.test <- X.full.test[,1:p0]
    Z1.test <- X.full.test[,(p0+1):(p0+p1)]
    Z.test <- cbind(Z0.test, Z1.test)

    ## Calculate effective beta
    beta_e = calculate_effective_beta(beta, Z1, idx_treat, interaction_design=interaction_design)
    beta_e.test = calculate_effective_beta(beta, Z1.test, idx_treat, interaction_design=interaction_design)

    ## Sample the response
    Y = sample_Y_linear(X.full, p1, beta_e, interaction_strength=0)

    #############################
    ## EXACT knockoff analysis ##
    #############################
    cat(sprintf("Applying EXACT method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.exact <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, cluster.method=cluster.method)
    discoveries.exact <- out.exact$discoveries
    ## Define hypotheses
    env.list <- out.exact$clusters
    hypotheses <- define_hypotheses(env.list, beta_e, Z, X)
    hypotheses.test <- define_hypotheses(env.list, beta_e.test, Z.test, X)
    ## Assess results
    res.exact <- evaluate_results_grouped(discoveries.exact, hypotheses) %>% mutate(Method="Exact", Test=FALSE, q=fdr.nominal, offset=offset)
    res.exact.test <- evaluate_results_grouped(discoveries.exact, hypotheses.test) %>% mutate(Method="Exact", Test=TRUE, q=fdr.nominal, offset=offset)
    cat(sprintf("Done with EXACT method.\n"))
    flush.console()

    #############################
    ## NAIVE knockoff analysis ##
    #############################
    cat(sprintf("Applying NAIVE method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.naive <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, naive=TRUE, cluster.method=cluster.method)
    discoveries.naive <- out.naive$discoveries
    ## Define hypotheses
    vars.split <- out.naive$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    hypotheses.test <- define_hypotheses(env.list, beta_e.test, Z.test, X)
    ## Assess results
    res.naive <- evaluate_results_grouped(discoveries.naive, hypotheses) %>% mutate(Method="Naive", Test=FALSE, q=fdr.nominal, offset=offset)
    res.naive.test <- evaluate_results_grouped(discoveries.naive, hypotheses.test) %>% mutate(Method="Naive", Test=TRUE, q=fdr.nominal, offset=offset)
    cat(sprintf("Done with NAIVE method.\n"))
    flush.console()

    #############################
    ## SPLIT knockoff analysis ##
    #############################
    cat(sprintf("Applying SPLIT method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.split <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, split=TRUE, cluster.method=cluster.method)
    discoveries.split <- out.split$discoveries
    ## Define hypotheses
    vars.split <- out.split$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    hypotheses.test <- define_hypotheses(env.list, beta_e.test, Z.test, X)
    ## Assess results
    res.split <- evaluate_results_grouped(discoveries.split, hypotheses) %>% mutate(Method="Split", Test=FALSE, q=fdr.nominal, offset=offset)
    res.split.test <- evaluate_results_grouped(discoveries.split, hypotheses.test) %>% mutate(Method="Split", Test=TRUE, q=fdr.nominal, offset=offset)
    cat(sprintf("Done with SPLIT method.\n"))
    flush.console()

    ###############################
    ## GLOBAL knockoff analysis ##
    ###############################
    cat(sprintf("Applying GLOBAL method...\n"))
    flush.console()
    ## Compute list of discoveries
    out.KFglobal <- analysis_knockoffs(Z, X, Y, q=fdr.nominal, offset=offset, KFglobal=TRUE)
    discoveries.KFglobal <- out.KFglobal$discoveries
    ## Define hypotheses
    vars.split <- out.KFglobal$splits
    hypotheses <- define_hypotheses(vars.split, beta_e, Z, X)
    hypotheses.test <- define_hypotheses(env.list, beta_e.test, Z.test, X)
    ## Assess results
    res.KFglobal <- evaluate_results_grouped(discoveries.KFglobal, hypotheses) %>% mutate(Method="KFglobal", Test=FALSE, q=fdr.nominal, offset=offset)
    res.KFglobal.test <- evaluate_results_grouped(discoveries.KFglobal, hypotheses.test) %>% mutate(Method="KFglobal", Test=TRUE, q=fdr.nominal, offset=offset)
    cat(sprintf("Done with GLOBAL method.\n"))
    flush.console()

    ## Combine results
    out <- rbind(res.exact, res.exact.test, res.naive, res.naive.test, res.split, res.split.test, res.KFglobal, res.KFglobal.test)

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
