#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

n.small = 1000        ## Number of samples to analyze
num.int=2
seed = 1              ## Random seed (for knockoffs)

n.small = as.integer(args[1])
num.int = as.integer(args[2])
seed = as.integer(args[3])

suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))
source("../../i-modelx/skf.R")
options("width"=200)

## Fixed parameters
seed.analysis <- 2022
set.seed(seed.analysis)

N <- 5
family = "binomial"
fdr.nominal=0.1

header <- tibble(seed=seed, fdr.nominal=fdr.nominal, num.int=num.int)

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
ifile.2 <- sprintf("../blood_donation_data/knockoffs/dt_donor_mar042018_simple_N5_seed%d.csv", seed)
data.raw.2 <- read_delim(ifile.2, delim="\t", col_types=cols(.default=col_double()))
data.raw <- cbind(data.raw.2, data.raw.1) %>% as_tibble()

#########################
## Data pre-processing ##
#########################
79203impute_na <- function(df) {
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
                   `age<25rh` = as.numeric(current_age < 25),
                   `edu<16` = as.numeric(edu_years < 16),
                   student = occupSTUD,
                   `recency<12` = as.numeric(recency_in_month < 12),
                   )
            
## Remove useless columns
data <- data %>% select("donate","donationvol",
                        "test_group", starts_with("X_"), starts_with("Xk_"),
                        "bloodtype", "weight", "male", "married", "age<25", "edu<16", "student", "resident",
                        "rh_neg", "recency<12" )
            
## Make the data set smaller (DEBUG)
if(n.small < nrow(data)) {
    idx.keep <- sample(1:nrow(data), pmin(n.small,nrow(data)))
    data <- data[idx.keep,]
}

X <- data %>% select(starts_with("X_")) %>% as.matrix()
Xk <- data %>% select(starts_with("Xk_")) %>% as.matrix()
X.treat <- cbind(X,Xk)
Y <- data$donate
family <- "binomial"
X.covar <- data.matrix(data[,14:ncol(data)])

## Print summary of the loaded data
cat(sprintf("Data set summary: %d observations, %d treatments, %d covariates.\n", 
            nrow(X), ncol(X), ncol(X.covar)))
cat(sprintf("Number of Y=0: %d; number of Y=1: %d.\n", sum(Y==0), sum(Y==1)))


#########################################
## Interaction lasso analysis (masked) ##
#########################################

data.tmp <- data %>% select(-donationvol, -test_group)

idx_treat.o <- 1+(1:N)
idx_treat.k <- 1+((N+1):(2*N))
idx_treat <- c(idx_treat.o, idx_treat.k)
idx_covar.noint <- (max(idx_treat)+1):(max(idx_treat)+2)
idx_covar.int <- (max(idx_covar.noint)+1):ncol(data.tmp)

Y <- as.integer(data.tmp$donate)-1
X <- data.matrix(data.tmp[,idx_treat.o])
Xk <- data.matrix(data.tmp[,idx_treat.k])
Z.noint <- data.matrix(data.tmp[,idx_covar.noint])
Z.int <- data.matrix(data.tmp[,idx_covar.int])

###############################
## Apply the knockoff filter ##
###############################

res <- list()

## Apply the Selective Knockoff Filter
cat("\n----------------------------------------\n")
cat(sprintf("Starting selective analysis...\n"))
cat("----------------------------------------\n")
res$skf <- skf_analysis(Y, X, Xk, Z.noint, Z.int, partition, family = "binomial", fdr.nominal=fdr.nominal, num.int=num.int,
                        naive=FALSE, split=FALSE, vanilla=FALSE, cross.prior=FALSE)

## Apply the Selective Knockoff Filter
cat("\n----------------------------------------\n")
cat(sprintf("Starting selective analysis (cp)...\n"))
cat("----------------------------------------\n")
res$skf.cp <- skf_analysis(Y, X, Xk, Z.noint, Z.int, partition, family = "binomial", fdr.nominal=fdr.nominal, num.int=num.int,
                           naive=FALSE, split=FALSE, vanilla=FALSE, cross.prior=TRUE)

## Apply the Vanilla Knockoff Filter
cat("\n----------------------------------------\n")
cat(sprintf("Starting vanilla analysis...\n"))
cat("----------------------------------------\n")
res$vanilla <- skf_analysis(Y, X, Xk, Z.noint, Z.int, partition, family = "binomial", fdr.nominal=fdr.nominal, num.int=num.int,
                            naive=FALSE, split=FALSE, vanilla=TRUE)

## Apply the naive Selective Knockoff Filter
cat("\n----------------------------------------\n")
cat(sprintf("Starting naive analysis...\n"))
cat("----------------------------------------\n")
res$naive <- skf_analysis(Y, X, Xk, Z.noint, Z.int, partition, family = "binomial", fdr.nominal=fdr.nominal, num.int=num.int,
                          naive=TRUE, split=FALSE, vanilla=FALSE)

## Apply the split Selective Knockoff Filter
cat("\n----------------------------------------\n")
cat(sprintf("Starting split analysis...\n"))
cat("----------------------------------------\n")
res$split <- skf_analysis(Y, X, Xk, Z.noint, Z.int, partition, family = "binomial", fdr.nominal=fdr.nominal, num.int=num.int,
                          naive=FALSE, split=TRUE, vanilla=FALSE)

summaries <- tibble()

for(method in names(res)) {
    discoveries <- res[[method]]$discoveries
    hypotheses <- res[[method]]$hypotheses
    ## Print the discoveries
    cat(sprintf("Method %s made %d discoveries:\n", method, nrow(discoveries)))
    if(nrow(discoveries)>0) {
        df.discoveries <- discoveries %>% left_join(hypotheses, by = c("Treatment", "Partition", "Label")) %>% print()
    } else {
        df.discoveries <- discoveries %>% mutate(Causal.prop=NA, Beta.sd=NA, n.sub=NA)
    }
    ## Evaluate the results
    df.summary <- evaluate_results(discoveries, hypotheses) %>% mutate(Method=method)
    summaries <- rbind(summaries, df.summary)
    ## Save results
    out.file <- sprintf("results/discoveries_real_n%d_nint%d_seed%d_%s.csv", n.small, num.int, seed, method)
    df.discoveries %>% bind_cols(header) %>% mutate(Method=method) %>% write_delim(out.file, delim=" ")
    cat(sprintf("List of discoveries written to %s\n", out.file))
    out.file <- sprintf("results/summary_real_n%d_nint%d_seed%d_%s.csv", n.small, num.int, seed, method)
    df.summary %>% bind_cols(header) %>% write_delim(out.file, delim=" ")
    cat(sprintf("Summary of discoveries written to %s\n", out.file))
}


cat(sprintf("Summaries for all methods:\n"))
summaries
