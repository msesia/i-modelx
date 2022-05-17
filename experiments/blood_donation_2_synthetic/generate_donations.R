#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

seed = 1
a = 1
b = a/2

a = as.double(args[1])
seed = as.integer(args[2])

suppressMessages(library(tidyverse))

###################
## Load the data ##
###################

mcoltypes <- cols(`test_group`=col_factor(),
                  `donate`=col_integer(),
                  `bloodtype`=col_factor(),
                  `rh_neg`=col_factor(),
                  `pass_test`=col_factor(),
                  `age_cat`=col_factor(),
                  `edu_cat`=col_factor(),
                  `occupation_cat`=col_factor(),
                  `test_group_new`=col_factor()
                 )

ifile.1 <- sprintf("../data/dt_donor_mar042018_simple.csv")
data.raw.1 <- read_delim(ifile.1, delim="\t",
                         col_types=mcoltypes)

ifile.2 <- sprintf("../data/knockoffs/dt_donor_mar042018_simple_N%d_seed%d.csv",5,seed)
data.raw.2 <- read_delim(ifile.2, delim="\t", col_types=cols())

data.raw <- cbind(data.raw.2, data.raw.1)

## Make the data set smaller (debug)
if(FALSE) {
    data.raw <- data.raw %>%
    slice(sample(n(), min(1000, n()))) %>%
    ungroup()
}

####################
## Pre-processing ##
####################

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
data <- data.raw %>%
    mutate(bloodtype = factor(bloodtype, c("O","A","B","AB")),
           rh_neg = factor(rh_neg, c(0,1)))

## Impute missing values
data <- impute_na(data)

## Make factors more meaningful
data <- data %>%
            mutate(donate = factor(donate, c(0,1)),
                   test_group = factor(test_group, c(0:6)),
                   age.below25 = as.numeric(current_age < 25),
                   education.below16 = as.numeric(edu_years < 16),
                   recency.below12 = as.numeric(recency_in_month < 12),
                   bloodtype.A = as.numeric(bloodtype=="A"),
                   bloodtype.B = as.numeric(bloodtype=="B"),
                   bloodtype.AB = as.numeric(bloodtype=="AB"),
                   )

## Remove useless columns
data <- data %>% select(-scramble_id, -test_group_new, -age_cat,
                        -eduyr16, -eduyr9, -eduyr16L, -eduyr12, -eduyrOther, -eduyr18, -eduyr12L, -edu_cat,
                        -last_donation_regtime, -pass_test, -formblood)

X <- data %>% select(starts_with("X_")) %>% as.matrix()
if(TRUE) {
    Y <- data$donate
    family <- "binomial"
} else {
    Y <- data$donationvol %>% as.numeric()
    family <- "gaussian"
}
Z <- data[,13:ncol(data)] %>% as.matrix()

cat(sprintf("Data set summary: %d observations, %d treatments, %d covariates.\n",
            nrow(X), ncol(X), ncol(Z)))
cat(sprintf("Number of Y=0: %d; number of Y=1: %d.\n", sum(Y==0), sum(Y==1)))

data.small <- data %>% select("donate", starts_with("X_"),
                              "male", "married", "resident", "age.below25",
                              "education.below16", "occupSTUD",
                              "rh_neg", "bloodtype.A", "bloodtype.B", "bloodtype.AB",
                              "recency.below12")

data.small %>% head()

###########################
## Define a causal model ##
###########################

compute_logistic_score <- function(x, a=0.4, b=0.4, c=0) {
    z <- -c +
         - 0.5*b * as.numeric(x["male"]) +
         + 0.5*b * as.numeric(x["married"]) +
         + 0.5*b * as.numeric(x["resident"]) +
         - 0.5*b * as.numeric(x["age.below25"]) +
         - 0.5*b * as.numeric(x["occupSTUD"]) +
         - 0.5*b * as.numeric(x["education.below16"]) +
         + 0.5*b * as.numeric(x["rh_neg"]) +
         - 0.5*b * (as.numeric(x["bloodtype.A"])+as.numeric(x["bloodtype.B"])+as.numeric(x["bloodtype.AB"])) +
         + a * as.numeric(x["X_1"] * (1-x["recency.below12"])) +
         + a * as.numeric(x["X_1"] * (1-x["resident"])) +
         + a * as.numeric(x["X_2"] * x["occupSTUD"]) +
         + a * as.numeric(x["X_2"] * x["occupSTUD"] * x["male"]) +
         + a * as.numeric(x["X_3"] * (1-x["male"])) +
         + a * as.numeric(x["X_3"] * (1-x["married"])) +
         + a * as.numeric(x["X_4"] * x["occupSTUD"]) +
         + a * as.numeric(x["X_4"] * x["male"] * x["occupSTUD"]) +
         + a * as.numeric(x["X_5"] * x["education.below16"])
    return(z)
}

compute_ite <- function(x, a, b, c) {
    ite <- rep(0, 5)
    ## T1
    ite[1] <- ite[1] + a * (1-x["recency.below12"])
    ite[1] <- ite[1] + a * (1-x["resident"])
    ## T2
    ite[2] <- ite[2] + a * (x["occupSTUD"])
    ite[2] <- ite[2] + a * (x["occupSTUD"]*x["male"])
    ## T3
    ite[3] <- ite[3] + a * (1-x["male"])
    ite[3] <- ite[3] + a * (1-x["married"])
    ## T4
    ite[4] <- ite[4] + a * (x["occupSTUD"])
    ite[4] <- ite[4] + a * (x["male"]*x["occupSTUD"])
    ## T5
    ite[5] <- ite[5] + a * (x["education.below16"])
    return(ite)
}

cat(sprintf("Estimating intercept for a=%.2f, b=%.2f... ", a, b))
idx.shuffle <- sample(1:nrow(data.small))
scores.pre <- sapply(1:1000, function(i.rel) {
    i <- idx.shuffle[i.rel]
    z <- compute_logistic_score(data.small[i,], a=a, b=b, c=0)
    return(z)
})
c = mean(scores.pre)
cat(sprintf("c=%.2f\n", c))

cat(sprintf("Computing logistic scores...\n"))
pb <- txtProgressBar(min = 0, max = nrow(data.small), style = 3)
scores <- sapply(1:nrow(data.small), function(i) {
    z <- compute_logistic_score(data.small[i,], a=a, b=b, c=c)
    setTxtProgressBar(pb, i)
    flush.console()
    return(z)
})
cat("\n")

cat(sprintf("Computing ITEs...\n"))
pb <- txtProgressBar(min = 0, max = nrow(data.small), style = 3)
ites <- t(sapply(1:nrow(data.small), function(i) {
    z <- compute_ite(data.small[i,], a=a, b=b, c=c)
    setTxtProgressBar(pb, i)
    flush.console()
    return(z)
}))
cat("\n")
df.ites <- tibble("X_1"=as.numeric(ites[,1]),
                  "X_2"=as.numeric(ites[,2]),
                  "X_3"=as.numeric(ites[,3]),
                  "X_4"=as.numeric(ites[,4]),
                  "X_5"=as.numeric(ites[,5]))

generate_Y <- function(scores) {
    invlogit = function(x) exp(x) / (1+exp(x))
    n <- length(scores)
    Y = rbinom(n, prob=invlogit(scores), size=1)
    return(Y)
}

set.seed(seed)
Y <- generate_Y(scores)

cat(sprintf("Donation rate: %.2f%%.",100*mean(Y)))
table(Y)

###################
## Save the data ##
###################

data.new <- data.raw %>% mutate(donate = Y)

scores.noisy <- scores + rnorm(length(scores), sd=0.2)
data.cont <- data.raw %>% mutate(donate = scores.noisy)

ofile <- sprintf("../data/synthetic/dt_donor_mar042018_synthetic_a%.1f_seed%d.csv", a, seed)
data.new %>% write_delim(ofile, delim="\t")
cat(sprintf("Synthetic data written to: %s\n", ofile))

ofile <- sprintf("../data/synthetic/dt_donor_mar042018_synthetic_continuous_a%.1f_seed%d.csv", a, seed)
data.cont %>% write_delim(ofile, delim="\t")
cat(sprintf("Synthetic data (continuous) written to: %s\n", ofile))

ofile <- sprintf("../data/synthetic/dt_donor_mar042018_synthetic_a%.1f_seed%d_ite.csv", a, seed)
df.ites %>% write_delim(ofile, delim="\t")
cat(sprintf("ITEs for synthetic data written to: %s\n", ofile))
