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

compute_interval <- function(X, Y, j, prob_treat, k, multivariate=FALSE, bound="lower", reference="auto",
                             n.steps = 1000, alpha=0.1, delta.max=20, step.size=1, verbose=FALSE) {

    Z = X[,j]
    n = length(Y)
    
    if(bound=="lower") {
        tail = "right"
        delta = -delta.max
    } else {
        tail = "left"
        delta = delta.max
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
        k = n-k
        ##delta = -delta
    }   

    if(multivariate) {
        ## Compute residuals
        mod = cv.glmnet(X[,-j],Y)
        Y = Y - predict(mod, X[,-j], s="lambda.min")
    }

    ## Debug What is the parameter of interest?
    ##delta_k = sort(beta[j,])[k]
    delta_k = NA
    if(verbose) cat(sprintf("%d-th (of %d) largest effect size: %.3f.\n", k, n, delta_k))

    history = cbind(rep(NA, n.steps+1), rep(NA, n.steps+1))
    history[1, ] <- c(0, delta)

    step = 1
    for(s in 1:n.steps) {
        step = step+1
        delta = interval.step(Y, Z, prob_treat, delta, k, reference=reference, tail=tail, multivariate=multivariate, alpha=alpha, step.size=step.size, verbose=verbose)
        history[s+1,] <- c(step, delta)
    }

    history = tibble(step=history[,1], bound=history[,2])

    ## Combine results
    out <- c()
    bound.data <- tail(history, round(nrow(history)/2))$bound
    bound.mean <- mean(bound.data)
    bound.sd <- sd(bound.data)    
    out$bound <- ifelse(bound=="lower", bound.mean-bound.sd, bound.mean+bound.sd)
    out$history <- history
    out$delta_k <- delta_k
    
    return(out)
}
