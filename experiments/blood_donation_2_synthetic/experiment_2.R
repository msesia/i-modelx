#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

n.small = 10000   ## How many samples to analyze?
a = 1.0          ## Interaction strength
seed = 1
model.class <- "glmnet"
ci.type <- "quantile"
tail <- "right"
k.rel <- 90
j <- 1

#####################
## Input arguments ##
#####################
if(TRUE) {
    n.small = as.integer(args[1])
    a = as.double(args[2])
    seed = as.integer(args[3])
    model.class = as.character(args[4])
    ci.type = as.character(args[5])
    tail <- as.character(args[6])
    k.rel = as.numeric(args[7])
    j = as.numeric(args[8])
}


suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))
options("width"=200)

env1 <- new.env()
sys.source("utils_ci_sharp.R", env1)

env2 <- new.env()
sys.source("utils_ci.R", env2)

## Fixed parameters
seed.analysis <- seed
set.seed(seed.analysis)

family = "gaussian"
reference <- "auto"
alpha=0.1

## Output file
out.dir = "results2"
out.file = sprintf("%s/n%d_a%1.f_%s_%s_t%s_k%d_seed%d_j%d.txt", out.dir, n.small, a, model.class, ci.type, tail, k.rel, seed, j)

## Set the quantile
k <- round(n.small*k.rel/100)

## Decide how many treatments to analyze
N = 5

tail.sign <- ifelse(tail=="right", "<=", ">=")
bound = ifelse(tail=="right", "lower", "upper")
header <- tibble(a=a, n=n.small, seed=seed,
                 alpha=alpha, k.rel=k.rel, k=k, tail.sign=tail.sign, model.class=model.class, ci.type=ci.type)

###################
## Load the data ##
###################

mcoltypes <- cols(`test_group`=col_factor(),
                  `donate`=col_double(),
                  `bloodtype`=col_factor(),
                  `pass_test`=col_factor(),
                  `age_cat`=col_factor(),
                  `edu_cat`=col_factor(),
                  `occupation_cat`=col_factor(),
                  `test_group_new`=col_factor()
                 )

ifile <- sprintf("../data/synthetic/dt_donor_mar042018_synthetic_continuous_a%.1f_seed%d.csv", a, seed)
data.raw <- read_delim(ifile, delim="\t", col_types=mcoltypes)

ifile.ite <- sprintf("../data/synthetic/dt_donor_mar042018_synthetic_a%.1f_seed%d_ite.csv", a, seed)
ite <- data.matrix(read_delim(ifile.ite, delim="\t", col_types=mcoltypes))


#########################
## Data pre-processing ##
#########################
impute_na <- function(df) {
   calc_mode <- function(x) {
      # List the distinct / unique values
      distinct_values <- unique(x)
      # Count the occurrence of each distinct value
      distinct_tabulate <- tabulate(match(x, distinct_values))
      # Return the value with the highest occurrence
      distinct_values[which.max(distinct_tabulate)]
    }

    df.out <- df
    for(cn in colnames(df)) {
        if(is.numeric(df.out[[cn]])) {
            df.out[[cn]][is.na(df.out[[cn]])] <- median(df.out[[cn]],na.rm=TRUE)
        }
        if(is.factor(df.out[[cn]])) {
            df.out[[cn]][is.na(df.out[[cn]])] <- calc_mode(df.out[[cn]])
        }
    }
    return(df.out)
}

## Remove empty strings
data.raw2 <- data.raw %>%
    mutate(bloodtype = factor(bloodtype, c("O","A","B","AB")))

## Impute missing values
data.raw2 <- impute_na(data.raw2)

## Make factors more meaningful
data <- data.raw2 %>% 
            mutate(#donate = factor(donate, c(0,1)),
                   test_group = factor(test_group, c(0:6)),
                   `age<25` = as.numeric(current_age < 25),
                   `edu<16` = as.numeric(edu_years < 16),
                   student = occupSTUD,
                   `recency<12` = as.numeric(recency_in_month < 12),
                   )
            
## Remove useless columns
data <- data %>% select("donate","donationvol",
                        "test_group", starts_with("X_"), starts_with("Xk_"),
                        "bloodtype", "weight", "male", "married", "age<25", "edu<16", "student", "resident",
                        "rh_neg", "recency<12" )

cat(sprintf("Loaded data set with %d rows.\n", nrow(data)))
    
## Make the data set smaller (DEBUG)
if(n.small < nrow(data)) {
    idx.keep <- sample(1:nrow(data), pmin(n.small,nrow(data)))
    data <- data[idx.keep,]
    ite <- ite[idx.keep,]
}

cat(sprintf("Sub-sampled data set to %d rows.\n", nrow(data)))

n <- nrow(data)
beta_e <- t(ite)
X <- data %>% select(starts_with("X_")) %>% as.matrix()
Y <- data$donate
X.covar <- data.matrix(data[,14:ncol(data)])

## Print summary of the loaded data
cat(sprintf("Data set summary: %d observations, %d treatments, %d covariates.\n", 
            nrow(X), ncol(X), ncol(X.covar)))
cat(sprintf("Mean of Y=0: %.2f; SD of Y=1: %.2f.\n", mean(Y), sd(Y)))

#########################################
## Interaction lasso analysis (masked) ##
#########################################

data.tmp <- data %>% select(-donationvol, -test_group, -starts_with("Xk_"))

idx_treat <- 1+(1:N)
idx_covar.noint <- (max(idx_treat)+1):(max(idx_treat)+2)
idx_covar.int <- (max(idx_covar.noint)+1):ncol(data.tmp)

Y <- data.tmp$donate
X <- data.matrix(data.tmp[,idx_treat])
Z.noint <- data.matrix(data.tmp[,idx_covar.noint])
Z.int <- data.matrix(data.tmp[,idx_covar.int])

####################################
## Construct CRT intervals        ##
####################################

conditional_prob_treat <- function(X, j) {
    n <- nrow(X)
    prob.treat <- rep(0,n)
    if(j==1) { # Reminder
        idx_0 <- which(rowSums(X[,-1])==0)
        prob.treat[idx_0] <- 11/25
        idx_1 <- which(rowSums(X[,-1])>0)
        prob.treat[idx_1] <- 1
    } else if(j==2) { # Individual reward
        prob.treat[which(apply(X[,-2], 1, function(x) return(all(x == c(1,0,0,0)))))] <- 0.5 # No other rewards
        prob.treat[which(apply(X[,-2], 1, function(x) return(all(x == c(1,1,0,1)))))] <- 1 # Friend + small
        prob.treat[which(apply(X[,-2], 1, function(x) return(all(x == c(1,1,1,0)))))] <- 0 # Friend + group reward
        prob.treat[which(apply(X[,-2], 1, function(x) return(all(x == c(1,1,0,0)))))] <- 0.5 # Friend
    } else if(j==3) { # Friend
        prob.treat[which(apply(X[,-3], 1, function(x) return(all(x == c(1,0,0,0)))))] <- 0.5 # No other rewards
        prob.treat[which(apply(X[,-3], 1, function(x) return(all(x == c(1,1,0,0)))))] <- 0.5 # Ind reward
        prob.treat[which(apply(X[,-3], 1, function(x) return(all(x == c(1,0,1,0)))))] <- 1 # Group rewards
        prob.treat[which(apply(X[,-3], 1, function(x) return(all(x == c(1,1,0,1)))))] <- 1 # Individual + small group
    } else if(j==4) { # Group reward
        prob.treat[which(apply(X[,-4], 1, function(x) return(all(x == c(1,0,1,0)))))] <- 0.5 # Friend
    } else if(j==5) { # Small group gift
        prob.treat[which(apply(X[,-5], 1, function(x) return(all(x == c(1,1,1,0)))))] <- 0.5 # Friend + reward
    }
    return(prob.treat)    
}

#####################
## DEBUG (dummies) ##
#####################

if(FALSE) {
    
    check_dummies <- function(X, j) {
        x <- X[,j]
        prob_treat <- conditional_prob_treat(X, j)
        xk <- rbinom(length(x), 1, prob_treat)
        cat("Marginal tables:\n")
        print(table(x))
        print(table(xk))
        X.full.1 <- X
        X.full.2 <- X.full.1
        X.full.2[,j] <- xk
        cat("Original covariance matrix:\n")
        print(cov(X.full.1))
        cat("Swapped covariance matrix:\n")
        print(cov(X.full.2))
    }

    check_dummies(X, 2)

}


#################
## Experiments ##
#################
results = tibble()

if(ci.type=="sharp") {
    multivariate.list <- c(FALSE)
    subset.list <- c(FALSE)
} else {
    multivariate.list <- c(FALSE, TRUE)
    subset.list <- c(FALSE)
}

for(multivariate in multivariate.list) {
    if(multivariate) {
        if(ci.type=="sharp") {
            interactions.list <- c(FALSE)
        } else {
            interactions.list <- c(FALSE, TRUE)
        }
    } else {
        interactions.list <- c(FALSE)
    }
    
    for(interactions in interactions.list) {
        for(subset in subset.list) {

            cat(sprintf("Running experiment (treatment variable %d, subset? %s, reference? %s, multivariate? %s, interactions? %s)...\n",
                        j, subset, reference, multivariate, interactions))

            ## Calculate conditional propensity scores
            prob_treat <- conditional_prob_treat(X, j)

            ## Combine the data
            X.full <- cbind(X, Z.noint, Z.int)

            ## Compute the sharp CI
            if(ci.type=="sharp") {
                if(model.class=="glmnet.fast") {
                    tmp <- env1$compute.crt.ci(X.full, Y, j, prob_treat, k, fast=TRUE, n.steps = 2500, step.size=0.5, alpha=alpha,
                                               model.class="glmnet", beta.max=5, beta.min=-5)
                } else {
                    tmp <- env1$compute.crt.ci(X.full, Y, j, prob_treat, k, fast=FALSE, n.steps = 2500, step.size=0.5, alpha=alpha,
                                               model.class=model.class, beta.max=5, beta.min=-5)
                }
                tmp.check <- env2$check_hypothesis(beta_e, j, round(0.5*n), 0, "lower", verbose=FALSE)
                ci.res <- c()
                ci.res$bound <- c("Lower", "Upper")
                ci.res$bound.value <- tmp$interval
                ci.res$tau_k <- tmp.check$tau_k
                ci.res$is_null <- all(beta_e[j,]==0)
            } else {
                tmp.check <- env2$check_hypothesis(beta_e, j, k, 0, tail, verbose=FALSE)
                tmp = env2$compute_interval(X.full, Y, j, prob_treat, k.rel, subset=subset, multivariate=multivariate, interactions=interactions,
                                            bound=bound, reference=reference, delta.max=5,
                                            n.steps = 5000, step.size=0.05, alpha=alpha, verbose=FALSE)
                ci.res <- c()
                ci.res$bound <- c(bound)
                ci.res$bound.value <- tmp$bound
                ci.res$tau_k <- tmp.check$tau_k
                ci.res$is_null <- tmp.check$is_null
            }

            res = tibble(treatment=j,
                         bound=ci.res$bound, bound.value=ci.res$bound.value,
                         tau_k=ci.res$tau_k, subset=subset, multivariate=multivariate, interactions=interactions,
                         reference=reference, null=ci.res$is_null)

            ## Store results
            results = rbind(results, res)
            cbind(results,header) %>%
                mutate_if(is.numeric, round, 4) %>%
                write_delim(out.file, delim=" ")

            print(results)

            cat(sprintf("Written results to: %s\n", out.file))
            flush.console()

        }
    }
}

cat(sprintf("All experiments completed.\n"))
flush.console()
