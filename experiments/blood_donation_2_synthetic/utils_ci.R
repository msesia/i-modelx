suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))

## Check hypothesis
check_hypothesis <- function(beta, j, k, delta, tail, verbose=FALSE) {
    # What is the parameter of interest?
    tau_k = sort(beta[j,])[k]
    if(verbose) {
        cat("Checking hypothesis:\n")
        cat(sprintf("%d-th (of %d) largest effect size: %.3f.\n", k, n, tau_k))
    }
    if(tail=="right") {
        tail.sign = "<="
        is_null = tau_k <= delta
    } else {
        tail.sign = ">="
        is_null = tau_k >= delta
    }
    if(verbose) {
        cat(sprintf("Null hypothesis: tau_(%d) %s %.3f\n", k, tail.sign, delta))
        cat(sprintf("Is the null hypothesis true? %s\n", is_null))
    }
    out = c()
    out$is_null = is_null
    out$tau_k = tau_k
    out$is_causal = any(beta[j,]!=0)
    return(out)
}

## Test statistics
phi_wilcoxon <- function(r) {
    return(r)
}

phi_stephenson <- function(r, s=10) {
    out = choose(r-1, s-1)
    out[r<s] = 0
    return(out)
}

t_stats <- function(z, y, side="right") {
    idx_treat = which(z==1)
    r = rank(y,ties.method="random")[idx_treat]
    if(side=="left") {
        r = length(y) - r
    }
    phi_vals = phi_wilcoxon(r)
    #phi_vals = phi_stephenson(r)
    T = log(sum(phi_vals))
    return(T)
}

## Transform Y to test one-sided hypotheses
adjust_Y <- function(Y, Z, Z_new, delta, k) {
    idx_treated = which(Z_new == 1)
    n = length(Z_new)
    m = length(idx_treated)
    l = min(m, n-k)
    I_k = idx_treated[which(rank(Y[idx_treated],ties.method="random")<=l)]
    eta = rep(delta,n)
    eta[I_k] = Inf
    delta_Y = (Z_new-Z)*eta
    ##delta_Y[Z_new==Z] = 0 ## DEBUG: make sure this is optimal
    delta_Y[which(is.nan(delta_Y))] = Inf
    Y_new = Y + delta_Y
    return(Y_new)
}

##############
## P-values ##
##############

compute_pval <- function(X, Y, j, prob_treat, k, delta, multivariate=FALSE, tail="right", reference="auto", K=100) {
    # Extract the treatment variable of interest
    Z = X[,j]
    n <- length(Y)

    if(multivariate) {
        ## Compute residuals
        mod = cv.glmnet(X[,-j],Y,alpha=0)
        Y = Y - predict(mod, X[,-j], s="lambda.min")
    } else {
        Y = Y
    }

    if (reference=="auto") {
        # Decide which reference to use to maximize power
        if(sum(prob_treat)<n/2) {
            reference = "control"
        } else {
            reference = "treated"
        }
    }

    if(reference=="control") {
        # Flip sign of response and treatment label
        Y = -Y
        Z = 1-Z
        prob_treat = 1 - prob_treat
    }

    # Flip response to test a left-tailed hypothesis
    if(tail=="left") {
        Y = -Y
        n = length(Y)
        k = n-k
        delta = -delta
        print("flip")
    }

    # Compute the null statistics
    T = t_stats(Z, Y)

    # Compute the test statistics
    T_null = sapply(1:K, function(b2) {
        Z_new = rbinom(n, 1, prob_treat)
        Y_new = adjust_Y(Y, Z, Z_new, delta, k)
        T_new = t_stats(Z_new, Y_new)
    })

    # Compute p-value
    pval = (1+sum(T_null>=T)) / (1+length(T_null))

    return(pval)
}

##########################
## Confidence intervals ##
##########################

interval.step <- function(Y, Z, prob_treat, delta, k, alpha=0.1, step.size=1, step=1,
                          multivariate = FALSE, verbose=FALSE, tail="right", reference="auto") {

    n <- length(Y)

    ## Compute statistics on original data
    T = t_stats(Z, Y)

    ## Sample a synthetic treatment variable
    Z_new = rbinom(n, 1, prob_treat)

    ## Computed synthetic Y
    if(tail=="left") {
        Y_new = adjust_Y(Y, Z, Z_new, -delta, k)
    } else {
        Y_new = adjust_Y(Y, Z, Z_new, delta, k)
    }

    ## Compute statistics on new data
    T_new = t_stats(Z_new, Y_new)

    ## Compute change in hypothesis
    if(tail=="right") {
        if(T_new < T) { ## Most common scenario (ideally happening with prob 1-alpha)
            ## Cautiously move up
            move <- step.size*(alpha)/(1+step/100)
        } else {
            ## Move back down
            move <- -step.size*(1-alpha)/(1+step/100)
        }
    } else {
         if(T_new < T) { ## Most common scenario (ideally happening with prob 1-alpha)
            ## Cautiously move up
            move <- -step.size*(alpha)/(1+step/100)
        } else {
            ## Move back down
            move <- step.size*(1-alpha)/(1+step/100)
        }
    }

    if(verbose) {
        cat(sprintf("Statistics: %.3f, %.3f. Delta = %.3f, Change = %.3f\n", T, T_new, delta, move))
    }

    delta <- delta + move
    return(delta)
}

predict_int_lasso_loo <- function(Y, X, X.covar, j, family="gaussian", verbose=FALSE) {
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
    str.treat.1 <- paste(paste("X",setdiff(1:N,j), sep="_"), collapse=" + ")
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
    pred <- predict(cv.fit.int, X.design, s="lambda.min")
    return(pred)
}

predict_lasso_loo <- function(Y, X, X.covar, j, family="gaussian", verbose=FALSE) {
    ## Create data frame with all data
    data.tmp <- data.frame(cbind(Y, X[,-j], X.covar))
    colnames(data.tmp)[1] <- "Y"
    idx_treat <- 1+(1:(N-1))
    idx_covar.noint <- (max(idx_treat)+1):(ncol(data.tmp))
    str.treat.1 <- paste(paste("X",setdiff(1:N,j), sep="_"), collapse=" + ")
    str.treat <- sprintf("%s", str.treat.1)
    str.noint <- paste(colnames(data.tmp)[idx_covar.noint], collapse=" + ")
    str.frmla <- sprintf("Y ~ %s + %s", str.treat, str.noint)
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
    pred <- predict(cv.fit.int, X.design, s="lambda.min")
    return(pred)
}

compute_interval <- function(X, Y, j, prob_treat, k.rel,
                             subset=FALSE, multivariate=FALSE, interactions=FALSE, bound="lower", reference="auto",
                             n.steps = 1000, delta.max=20, alpha=0.1, step.size=1, verbose=FALSE) {

    n <- length(Y)

    if(subset) {
        idx_const <- which((prob_treat==0)|(prob_treat==1))
        idx_test <- setdiff(1:n, idx_const)
        idx_train <- 1:n
    } else {
        idx_train <- 1:n
        idx_test <- 1:n
    }

    n <- length(idx_test)
    k <- round(k.rel*n/100)

    Z = X[,j]

    if(bound=="lower") {
        tail = "right"
        delta = -delta.max
    } else {
        tail = "left"
        delta = delta.max
    }

    if (reference=="auto") {
        # Decide which reference to use to maximize power
        if(sum(prob_treat[idx_test])<n/2) {
            reference = "control"
        } else {
            reference = "treated"
        }
    }

    if(reference=="control") {
        # Flip sign of response and treatment label
        Y = -Y
        Z = 1-Z
        prob_treat = 1 - prob_treat
    }

    # Flip response to test a left-tailed hypothesis
    if(tail=="left") {
        Y = -Y
        k = n-k
    }
    if(multivariate) {
        ## Compute residuals
        X.treat = X[,1:5]
        X.covar = X[,6:ncol(X)]
        ##X.treat[idx_test,j] <- rbinom(length(idx_test), 1, prob_treat[idx_test])
        if(interactions) {
            pred.test = as.numeric(predict_int_lasso_loo(Y, X.treat, X.covar, j))
        } else {
            pred.test = as.numeric(predict_lasso_loo(Y, X.treat, X.covar, j))
        }
        Y <- Y - pred.test
    }

    ## Debug What is the parameter of interest?
    ##delta_k = sort(beta[j,])[k]
    delta_k = NA
    ##if(verbose) cat(sprintf("%d-th (of %d) largest effect size: %.3f.\n", k, n, delta_k))

    history = cbind(rep(NA, n.steps+1), rep(NA, n.steps+1))
    history[1, ] <- c(0, delta)

    step = 1
    for(s in 1:n.steps) {
        step = step+1
        delta = interval.step(Y[idx_test], Z[idx_test], prob_treat[idx_test],
                              delta, k, reference=reference, tail=tail, multivariate=multivariate, alpha=alpha, step.size=step.size, verbose=verbose)
        history[s+1,] <- c(step, delta)
    }

    history = tibble(step=history[,1], bound=history[,2])

    ## Combine results
    out <- c()
    bound.data <- tail(history, round(0.5*nrow(history)))$bound
    bound.mean <- mean(bound.data)
    bound.sd <- sd(bound.data)
    out$bound <- ifelse(bound=="lower", bound.mean-bound.sd, bound.mean+bound.sd)
    out$history <- history
    out$delta_k <- delta_k

    return(out)
}
