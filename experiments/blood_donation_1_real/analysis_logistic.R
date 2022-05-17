#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

n.small = 80000        ## Number of samples to analyze
seed = 1              ## Random seed (for knockoffs)

n.small = as.integer(args[1])
seed = as.integer(args[2])

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
ifile.2 <- sprintf("../blood_donation_data/knockoffs/dt_donor_mar042018_simple_N5_seed%d.csv", seed)
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

##############################################
## Logistic regression without interactions ##
##############################################

## With covariates (does not converge)
glm.fit <- glm(Y~., data, family="binomial")
summary(glm.fit)

## No covariates
glm.fit <- glm(Y~X_1+X_2+X_3+X_4+X_5, data, family="binomial")
summary(glm.fit)

pvals <- coef(summary(glm.fit))[-1,4]
pvals.adj <- p.adjust(pvals, method="BH")


#######################
## Marginal analysis ##
#######################

pvals.fisher <- sapply(1:N, function(j) fisher.test(X[,j], Y)$p.value)
round(pvals.fisher,3)
round(p.adjust(pvals.fisher, method="BH"),3)
