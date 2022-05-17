#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

n.small = 1000        ## Number of samples to analyze
j = 1                 ## Treatment to analyze

n.small = as.integer(args[1])
j = as.integer(args[2])

suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))
options("width"=200)

## Fixed parameters
seed.analysis <- 2022
set.seed(seed.analysis)

N <- 5
family = "binomial"
fdr.nominal=0.1

###################
## Load the data ##
###################

mcoltypes <- cols(`test_group`=col_factor(),
                  `donate`=col_integer(),
                  `bloodtype`=col_factor(),
                  `pass_test`=col_factor(),
                  `age_cat`=col_factor(),
                  `edu_cat`=col_factor(),
                  `occupation_cat`=col_factor(),
                  `test_group_new`=col_factor()
                 )

ifile.1 <- sprintf("../blood_donation_data/dt_donor_mar042018_simple.csv")
data.raw.1 <- read_delim(ifile.1, delim="\t", col_types=mcoltypes)
ifile.2 <- sprintf("../blood_donation_data/knockoffs/dt_donor_mar042018_simple_N5_seed%d.csv", 1)
data.raw.2 <- read_delim(ifile.2, delim="\t", col_types=cols(.default=col_double()))
data.raw <- cbind(data.raw.2, data.raw.1) %>% as_tibble()

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
            mutate(donate = factor(donate, c(0,1)),
                   test_group = factor(test_group, c(0:6)),
                   `age<25` = as.numeric(current_age < 25),
                   `edu<16` = as.numeric(edu_years < 16),
                   student = occupSTUD,
                   `recency<12` = as.numeric(recency_in_month < 12),
                   )
            
## Remove useless columns
data <- data %>% select("donate","donationvol",
                        "test_group", starts_with("X_"),
                        "bloodtype", "weight", "male", "married", "age<25", "edu<16", "student", "resident",
                        "rh_neg", "recency<12" )
            
## Make the data set smaller (DEBUG)
if(n.small < nrow(data)) {
    idx.keep <- sample(1:nrow(data), pmin(n.small,nrow(data)))
    data <- data[idx.keep,]
}

X <- data %>% select(starts_with("X_")) %>% as.matrix()
Y <- data$donate
family <- "binomial"
X.covar <- data.matrix(data[,9:ncol(data)])

## Print summary of the loaded data
cat(sprintf("Data set summary: %d observations, %d treatments, %d covariates.\n", 
            nrow(X), ncol(X), ncol(X.covar)))
cat(sprintf("Number of Y=0: %d; number of Y=1: %d.\n", sum(Y==0), sum(Y==1)))

###################################
## Conditional propensity scores ##
###################################

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


predict_int_lasso <- function(Y, X, X.covar, j, verbose=FALSE) {
    ## Create data frame with all data
    data.tmp <- data.frame(cbind(Y, X[,-j], X.covar))
    colnames(data.tmp)[1] <- "Y"
    idx_treat <- 1+(1:(N-1))
    idx_covar.noint <- (max(idx_treat)+1):(max(idx_treat)+2)
    idx_covar.int <- (max(idx_covar.noint)+1):ncol(data.tmp)
    ## Make list of interactions
    interactionCandidates <- idx_treat
    interactionPairs <- lapply(idx_treat, function(i) {return(cbind(i,idx_covar.int))})
    interactionPairs <- as_tibble(do.call("rbind", interactionPairs))
    df.int <- interactionPairs %>% mutate(Treatment = colnames(data.tmp)[i],
                                          Covariate = paste("`",colnames(data.tmp)[idx_covar.int],"`",sep="")) %>%
        select(-i, -idx_covar.int)
    ## Define the formula
    str.treat.1 <- paste(paste("X",setdiff(1:N,j),sep="_"), collapse=" + ")
    str.treat <- sprintf("%s", str.treat.1)
    str.noint <- paste(colnames(data.tmp)[idx_covar.noint], collapse=" + ")
    str.int <- paste(df.int$Treatment, "*", df.int$Covariate, collapse=" + ")
    str.frmla <- sprintf("Y ~ %s + %s + %s", str.treat, str.noint, str.int)
    if(verbose) {
        cat(sprintf("Lasso formula:\n"))
        print(str.frmla)
    }
    ## Create the data matrix
    frmla <- as.formula(str.frmla)
    X.design <- model.matrix(frmla, data.tmp, keep.order = TRUE)[,-1]
    ## Fit the lasso with interactions
    cv.fit.int <- cv.glmnet(X.design, Y, family=family, alpha=1)
    beta.hat.int <- coef(cv.fit.int, s="lambda.min")
    ## Make predictions
    pred <- predict(cv.fit.int, X.design)
    return(pred)
}

compute_pval_sharp <- function(Y, X, X.covar, j, family="gaussian", n_mc=100, verbose=TRUE) {
    # Extract the treatment variable of interest
    Z = X[,j]
    ## Calculate conditional propensity scores
    prob_treat <- conditional_prob_treat(X, j)   
    ## Interaction lasso analysis (LOO)
    if(verbose) {
        cat("Fitting lasso model... ")
    }
    if(verbose) {
        cat("Done.\n")
    }
    pred <- predict_int_lasso(Y, X, X.covar, j, verbose=FALSE)
    t_stats <- function(Y, x, pred) {
        mod.fit <- glm(Y~pred+x, family=family)
        as.numeric(abs(summary(mod.fit)$coef[,3][2]))
    }    
    # Compute the null statistics
    x <- X[,j]
    T = t_stats(Y, x, pred)
    # Compute the test statistics
    if(verbose) {
        cat(sprintf("Computing CRT with %d randomizations... ", n_mc))
    }
    pb <- txtProgressBar(min = 0, max = n_mc, style = 3)
    T_null = sapply(1:n_mc, function(b2) {
        x_new = rbinom(length(prob_treat), 1, prob_treat)
        T_new = t_stats(Y, x_new, pred)
        setTxtProgressBar(pb, b2)
        flush.console()
        return(T_new)
    })
    if(verbose) {
        cat("Done.\n")
    }
    # Compute p-value
    pval = (1+sum(T_null>=T)) / (1+length(T_null))
    return(pval)
}

results <- tibble()
cat(sprintf("Computing CRT p-value for treatment %d...\n", j))
pval <- compute_pval_sharp(Y, X, X.covar, j, family=family, n_mc=10000)
cat(sprintf("P-value = %.3g\n", pval))
res <- tibble(n=n.small, treatment=j, pval=pval)
results <- rbind(results, res)
out.file <- sprintf("results_crt/pvalues_real_n%d_j%d.csv", n.small, j)
cat(sprintf("List of p-values written to %s\n", out.file))
results %>% write_delim(out.file, delim=" ")
