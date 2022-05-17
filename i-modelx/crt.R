suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))

## Check hypothesis
check_hypothesis = function(beta, j, k, delta, tail, verbose=FALSE) {
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

t_stats = function(X, j, Z, Y) {
    X[,j] = Z
    mod = cv.glmnet(X,Y,alpha=1)
    T = abs(coef(mod, s="lambda.min")[-1])[j]
    return(T)
}

##############
## P-values ##
##############

compute_pval_sharp = function(X, Y, j, prob_treat, delta) {
    # Extract the treatment variable of interest
    Z = X[,j]

    # Shift response
    Y = Y - delta * Z

    # Compute the null statistics
    T = t_stats(X, j, Z, Y)

    # Compute the test statistics
    T_null = sapply(1:100, function(b2) {
        Z_new = rbinom(n, 1, prob_treat)
        T_new = t_stats(X, j, Z_new, Y)
    })

    # Compute p-value
    pval = (1+sum(T_null>=T)) / (1+length(T_null))

    return(pval)
}
