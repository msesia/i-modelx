suppressMessages(library(tidyverse))
suppressMessages(library(randomForest))

fit_int_lasso <- function(Y, X, X.covar, j, family="gaussian", verbose=FALSE) {
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
    pred <- predict(cv.fit.int, X.design)
    return(pred)
}

check.convergence <- function(history, alpha=0.1, roll.width=250, beta.min=-10, beta.max=10) {
    n.steps <- nrow(history) - 1
    ## Convergence diagnostics
    alpha.se <- sqrt(alpha*(1-alpha)) / sqrt(roll.width)
    history.delta <- diff(history)
    diagnostics <- tibble(Step=seq(1,n.steps),
                          Lower=history[-1,2], Upper=history[-1,3],
                          Delta.lower=history.delta[,2], Delta.upper=history.delta[,3]) %>%
        mutate(Change.lower = sign(Delta.lower)<0, Change.upper = sign(Delta.upper)>0) %>%
        mutate(Lower.mean = zoo::rollmean(Lower, roll.width, na.pad=TRUE),
               Upper.mean = zoo::rollmean(Upper, roll.width, na.pad=TRUE),
               RL = zoo::rollmean(Change.lower, roll.width, na.pad=TRUE),
               RU = zoo::rollmean(Change.upper, roll.width, na.pad=TRUE)) %>%
        drop_na(RL, RU) %>%
        select(-Change.lower, -Change.upper, -Delta.lower, -Delta.upper) %>%
        mutate(Converged.A = (abs(RL-alpha)<2*alpha.se)*(abs(RU-alpha)<2*alpha.se)==1,
               Converged.B = abs(Lower.mean-Upper.mean)<0.01,
               Converged = any(Converged.A, Converged.B),
               Unfeasible.L = Lower < beta.min,
               Unfeasible.U = Upper > beta.max,
               Failed = (RL>alpha+2*alpha.se) && (RU>alpha+2*alpha.se),
               Converged = ifelse(Unfeasible.L|Unfeasible.U, FALSE, Converged),
               Failed = ifelse(Unfeasible.L|Unfeasible.U, FALSE, Failed))

    diagnostics <- tail(diagnostics, n=1)
    if(diagnostics$Converged) {
        status <- "converged"
    } else {
        if(diagnostics$Failed) {
            status <- "failed"
        } else if(diagnostics$Unfeasible.L) {
            if(diagnostics$Unfeasible.U) {
                status <- "unfeasible.both"
            } else {
                status <- "unfeasible.lower"
            }
        } else if(diagnostics$Unfeasible.U) {
            status <- "unfeasible.upper"
        } else {
            status <- "incomplete"
        }
    }
    return(status)
}

plot.history <- function(history, beta0.true=0, alpha=0.1, roll.width=250, beta.min=-10, beta.max=10) {
    n.steps <- nrow(history) - 1

    p1 <- tibble(Step=seq(0,n.steps), Lower=history[,2], Upper=history[,3]) %>%
        gather(Lower, Upper, key="Limit", value="Value") %>%
        ggplot(aes(x=Step, y=Value, color=Limit, group=Limit)) +
        geom_point() +
        geom_hline(yintercept=beta0.true, linetype=2) +
        theme_bw()

    alpha.se <- sqrt(alpha*(1-alpha)) / sqrt(roll.width)
    history.delta <- diff(history)
    p2 <- tibble(Step=seq(1,n.steps), Lower=history.delta[,2], Upper=history.delta[,3]) %>%
        mutate(Change.lower = sign(Lower)<0, Change.upper = sign(Upper)>0) %>%
        mutate(Rate.lower = zoo::rollmean(Change.lower, roll.width, na.pad=TRUE),
               Rate.upper = zoo::rollmean(Change.upper, roll.width, na.pad=TRUE)) %>%
        gather(Rate.lower, Rate.upper, key="Limit", value="Rate") %>%
        ggplot(aes(x=Step, y=Rate, color=Limit, group=Limit)) +
        geom_line() +
        geom_hline(yintercept=alpha, linetype=2) +
        geom_hline(yintercept=alpha-2*alpha.se, linetype=3) +
        geom_hline(yintercept=alpha+2*alpha.se, linetype=3) +
        ylim(0,1) +
        theme_bw()

    p3 <- tibble(Step=seq(0,n.steps), Lower=history[,2], Upper=history[,3]) %>%
        mutate(Lower = zoo::rollmean(Lower, roll.width, na.pad=TRUE),
               Upper = zoo::rollmean(Upper, roll.width, na.pad=TRUE)) %>%
        gather(Lower, Upper, key="Limit", value="Value") %>%
        ggplot(aes(x=Step, y=Value, color=Limit, group=Limit)) +
        geom_point() +
        geom_hline(yintercept=beta0.true, linetype=2) +
        ylim(beta.min,beta.max) +
        theme_bw()

    gridExtra::grid.arrange(p1, p2, p3, ncol=1)
}

crt.loop <- function(j, X, prob_treat, y, beta.lower, beta.upper, beta0.grid, lamb.precomputed,
                     n.steps=250, roll.width=250, step.size=1, step.start=0, alpha=0.1, fast=FALSE,
                     beta.min=-10, beta.max=10, model.class="glmnet") {

    is_inactive <- function(prob_treat) {
        w = prob_treat*(1-prob_treat)
        return(w < quantile(w,0.1))
    }

    n <- length(prob_treat)
    p <- ncol(X)

    glmnet.tune <- function(beta0, beta0.grid, lamb.precomputed) {
        beta.bin <- .bincode(beta0, beta0.grid)
        if(is.na(beta.bin)) beta.bin <- 1
        X.tmp <- X
        if(is.na(lamb.precomputed[beta.bin])) {
            ##cat(sprintf("Computing optimal lambda for bin %d...\n", beta.bin))
            z <- X[,j]
            y.bar <- y - z * beta0
            ## Pre-compute optimal regularization
            z.new <- rbinom(length(prob_treat), 1, prob_treat)
            X.tmp[,j] <- z.new
            y.tmp <- y.bar + z.new * beta0
            cv.lasso <- glmnet::cv.glmnet(X.tmp, y.tmp, family="gaussian", nlambda=50)
            lamb <- cv.lasso$lambda.min
            lamb.precomputed[beta.bin] <- lamb
        } else {
            ##cat(sprintf("Found optimal lambda for bin %d.\n", beta.bin))
            lamb <- lamb.precomputed[beta.bin]
        }
        res <- c()
        res$lamb <- lamb
        res$lamb.precomputed <- lamb.precomputed
        return(res)
    }

    compute.importance.glmnet <- function(X_tmp, y.bar, z, beta0, lamb=NULL, base.model=NULL) {
        X_tmp[,j] <- z
        y.tmp <- y.bar
        if(is.null(base.model)) {
            ## Penalty factor makes sure beta.hat is not zero
            penalty.factor <- rep(1, p)
            penalty.factor[j] <- 0
            ## Compute lasso
            if(lamb>1) {
                lamb.seq <- seq(lamb*10, lamb, length.out=20)
            } else {
                lamb.seq <- seq(1, lamb, length.out=20)
            }
            opt.lasso <- glmnet::glmnet(X_tmp, y.tmp, family="gaussian", lambda=lamb.seq, penalty.factor=penalty.factor)
            beta.hat <- coef(opt.lasso, s=lamb)[-1][j]
        } else {
            idx.active <- which(!is_inactive(prob_treat))
            pred <- as.numeric(X_tmp[idx.active,] %*% base.model)
            residuals <- y.tmp[idx.active]-pred
            beta.hat <- as.numeric(coef(lm(residuals~z[idx.active]))[2])
            beta.hat <- ifelse(is.na(beta.hat), 0, beta.hat)
        }
        imp.stat <- abs(beta.hat) + rnorm(1, sd=1e-6)
        return(imp.stat)
    }

    compute.importance.glmnet.inter <- function(X_tmp, y.bar, z, beta0, lamb=NULL, base.model=NULL) {                    
        X_tmp[,j] <- z
        pred <- as.numeric(fit_int_lasso(y.bar, X_tmp[,1:5], X_tmp[,6:ncol(X_tmp)], j, family="gaussian", verbose=FALSE))
        beta.hat <- as.numeric(coef(lm(y.bar~pred+z))[3])       
        beta.hat <- ifelse(is.na(beta.hat), 0, beta.hat)
        imp.stat <- abs(beta.hat) + rnorm(1, sd=1e-6)
        return(imp.stat)
    }

    compute.importance.rf <- function(X_tmp, y.bar, z, beta0, lamb=NULL, base.model=NULL) {
        X_tmp[,j] <- z
        y.tmp <- y.bar
        df = as_tibble(X_tmp) %>% mutate(Y=y.tmp)
        tree = ranger::ranger(Y~., df, importance = "impurity", num.trees=100, max.depth=2)
        vi = as.numeric(tree$variable.importance)
        imp.stat <- vi[j] + rnorm(1, sd=1e-6)
        return(imp.stat)
    }

    stats.continuous <- function(beta0, lamb=NULL, base.model=NULL) {
        ## Compute baseline response
        z <- X[,j]
        y.bar <- y - z * beta0

        if(model.class=="glmnet") {
            compute.importance = compute.importance.glmnet
        } else if(model.class=="glmnet.inter") {
            compute.importance = compute.importance.glmnet.inter
        } else if(model.class=="rf") {
            compute.importance = compute.importance.rf
        } else {
            stop("Error: unknown model class!")
        }

        ## Compute statistics on original data
        stat.original <- compute.importance(X, y.bar, z, beta0, lamb=lamb, base.model=base.model)

        ## Resample the treatment allocation
        z.new <- rbinom(length(prob_treat), 1, prob_treat)

        ## Compute statistics on perturbed data
        stat.new <- compute.importance(X, y.bar, z.new, beta0, lamb=lamb, base.model=base.model)

        ## Return results
        return(c(stat.original, stat.new))
    }

    crt.step <- function(beta0, alpha=0.1, step.size=1, upper=TRUE, step=1, lamb=NULL, base.model=NULL) {
        cat(sprintf("Step %d: beta0 = %.3f: ", step, beta0))
        ## Compute statistics on original and perturbed data
        T <- stats.continuous(beta0, lamb=lamb, base.model=base.model)
        t.ori <- T[1]
        t.new <- T[2]
        ## Compute change in hypothesis
        if(upper) {
            if(t.ori > t.new) { ## Most common scenario (ideally happening with prob 1-alpha)
                ## Cautiously move down
                delta <- -step.size*alpha / (1+step/100)
            } else {
                ## Move back up
                delta <- step.size*(1-alpha) / (1+step/100)
            }
        } else {
            if(t.ori > t.new) { ## Most common scenario (ideally happening with prob 1-alpha)
                ## Cautiously move up
                delta <- step.size*(alpha)/(1+step/100)
            } else {
                ## Move back down
                delta <- -step.size*(1-alpha)/(1+step/100)
            }
        }
        if(upper) {
            cat(sprintf("Upper: %.3f, %.3f. Delta = %.3f\n", T[1], T[2], delta))
        } else {
            cat(sprintf("Lower: %.3f, %.3f. Delta = %.3f\n", T[1], T[2], delta))
        }
        new.beta0 <- beta0 + delta
        return(new.beta0)
    }

    if(fast) {
        if(model.class=="glmnet") {
            idx.inactive <- which(is_inactive(prob_treat))
            if(length(idx.inactive)>2) {
                cv.lasso <- glmnet::cv.glmnet(X[idx.inactive,], y[idx.inactive], family="gaussian")
                lamb <- cv.lasso$lambda.min
                if(lamb>1) {
                    lamb.seq <- seq(lamb*10, lamb, length.out=20)
                } else {
                    lamb.seq <- seq(1, lamb, length.out=20)
                }
                base.lasso <- glmnet::glmnet(X[idx.inactive,], y[idx.inactive], family="gaussian", lambda=lamb.seq)
                base.model <- coef(base.lasso, s=lamb)[-1]
            } else {
                base.model <- rep(0, ncol(X))
            }
        } else if(model.class=="glmnet.inter") {
            idx.inactive <- which(is_inactive(prob_treat))
            data <- as.data.frame(cbind(y, X))
            colnames(data) <- c("Y", paste("Z", seq(1,p0), sep=""), paste("ZI", seq(1,p1), sep=""), paste("T", seq(1,p1), sep=""))
            f.str <- sprintf("Y ~ . + (%s):(%s)", paste(paste("ZI", seq(1,p1), sep=""), collapse="+"), paste(paste("T", seq(1,p1), sep=""), collapse="+"))
            f <- as.formula(f.str)
            x <- model.matrix(f, data)[, -1]
            x.names <- colnames(x)
            y.tmp <- data$Y
            cv.lasso <- glmnet::cv.glmnet(x, y.tmp, family="gaussian")
            lamb <- cv.lasso$lambda.min
            if(lamb>1) {
                lamb.seq <- seq(lamb*10, lamb, length.out=20)
            } else {
                lamb.seq <- seq(1, lamb, length.out=20)
            }
            base.lasso <- glmnet::glmnet(x[idx.inactive,], y.tmp[idx.inactive], family="gaussian", lambda=lamb.seq)
            base.model <- coef(base.lasso, s=lamb)[-1]
        } else {
            stop("Error: unknown model class!")
        }
    } else {
        base.model <- NULL
    }

    beta.seq <- cbind(rep(NA, n.steps), rep(NA, n.steps),rep(NA, n.steps))
    step.size.lower <- step.size
    step.size.upper <- step.size
    for(s in 1:n.steps) {
        step <- step.start + s
        ## Update lower bound
        if(fast) {
            beta.lower <- crt.step(beta.lower, alpha=alpha, step.size=step.size.lower, upper=FALSE, step=step,
                                   base.model=base.model)
        } else {
            cv.res.lower <- glmnet.tune(beta.lower, beta0.grid, lamb.precomputed)
            lamb.lower <- cv.res.lower$lamb
            lamb.precomputed <- cv.res.lower$lamb.precomputed
            beta.lower <- crt.step(beta.lower, alpha=alpha, step.size=step.size.lower, upper=FALSE, step=step,
                                   lamb=lamb.lower)
        }
        ## Update upper bound
        if(fast) {
            beta.upper <- crt.step(beta.upper, alpha=alpha, step.size=step.size.upper, upper=TRUE, step=step,
                                   base.model=base.model)
        } else {
            cv.res.upper <- glmnet.tune(beta.upper, beta0.grid, lamb.precomputed)
            lamb.upper <- cv.res.upper$lamb
            lamb.precomputed <- cv.res.upper$lamb.precomputed
            beta.upper <- crt.step(beta.upper, alpha=alpha, step.size=step.size.upper, upper=TRUE, step=step,
                                   lamb=lamb.upper)
        }
        ## Make sure lower < upper
        if(beta.lower > beta.upper) {
            tmp <- beta.upper
            beta.upper <- beta.lower
            beta.lower <- tmp
        }
        ## Make sure neither bound is too extreme
        beta.lower <- pmax(beta.min*1.25, beta.lower)
        beta.upper <- pmin(beta.max*1.25, beta.upper)
        beta.seq[s,] <- c(step, beta.lower, beta.upper)
    }
    return(beta.seq)
}

compute.crt.ci <- function(X, y, j, prob_treat, seed=seed, alpha=0.1, fast=FALSE,
                           step.size=1, n.steps=1000, roll.width=250, beta.min=-10, beta.max=10,
                           model.class="glmnet") {

    ## Check whether there is information in the data
    if(sum(prob_treat*(1-prob_treat))==0) {
        out <- c()
        out$interval <- c(beta.min, beta.max)
        out$estimate <- mean(out$interval)
        out$history <- tibble()
        out$status <- "SKIP"
        return(out)
    }

    set.seed(seed)

    beta0.grid <- seq(beta.min,beta.max,length.out=20)
    lamb.precomputed <- rep(NA, length(beta0.grid))

    running <- TRUE
    step.start <- 0
    beta.lower <- beta.min
    beta.upper <- beta.max
    history <- matrix(c(0,beta.lower,beta.upper),1,3)
    num.fail <- 0

    while(running) {
        ## Estimate lower and upper bounds
        beta.seq <- crt.loop(j, X, prob_treat, y, beta.lower, beta.upper, beta0.grid, alpha=alpha,
                             lamb.precomputed, n.steps=n.steps, roll.width=roll.width, fast=fast,
                             step.size=step.size, step.start=step.start, beta.min=beta.min, beta.max=beta.max,
                             model.class=model.class)
        history <- rbind(history, beta.seq)

        ## Check convergence
        status <- check.convergence(history, roll.width=roll.width, beta.min=beta.min, beta.max=beta.max)

        if(status=="converged") {
            ## Exit loop
            running <- FALSE
            cat(sprintf("Status check: CONVERGED!\n"))
        } else if(status=="incomplete") {
            ## Continue for another 'n.steps' steps
            step.start <- step.start + n.steps
            beta.lower <- tail(history[,2], n=1)
            beta.upper <- tail(history[,3], n=1)
            cat(sprintf("Status check: INCOMPLETE. Continuing...\n"))
        } else if (status=="unfeasible.both") {
            running <- FALSE
            cat(sprintf("Status check: CONVERGED (both bounds are unfeasible)!\n"))
        } else if (status=="unfeasible.lower") {
            running <- FALSE
            cat(sprintf("Status check: CONVERGED (lower bound is unfeasible)!\n"))
        } else if (status=="unfeasible.upper") {
            running <- FALSE
            cat(sprintf("Status check: CONVERGED (upper bound is unfeasible)!\n"))
        } else {
            num.fail <- num.fail + 1
            if(num.fail < 2) {
                cat(sprintf("Status check: FAILED. Starting over...\n"))
                step.size <- step.size/2
                step.start <- 0
            } else {
                cat(sprintf("Status check: FAILED too many times. Giving up...\n"))
                running <- FALSE
            }
            beta.lower <- beta.min
            beta.upper <- beta.max
            history <- matrix(c(0,beta.lower,beta.upper),1,3)
        }

        if(nrow(history)>10000) {
            stop("Error: time out!")
        }
    }

    ## Convergence diagnostics
    if(FALSE) {
        plot.history(history, roll.width=100, beta0.true=beta[j])
    }

    ## Return results
    out <- c()
    extract.interval <- function(history, status, delta=0, roll.width=250, alpha=0.1) {
        if(status=="failed") {
            lower <- history[1,2]
            upper <- history[1,3]
        } else {
            n.steps <- nrow(history) - 1
            history.tail <- tibble(Step=seq(0,n.steps), Lower=history[,2], Upper=history[,3]) %>%
                tail(roll.width)
            lower <- mean(history.tail$Lower) - delta
            upper <- mean(history.tail$Upper) + delta
            if (status=="unfeasible.both") {
                lower <- history[1,2]
                upper <- history[1,3]
            } else if (status=="unfeasible.lower") {
                lower <- history[1,2]
            } else if (status=="unfeasible.upper") {
                upper <- history[1,3]
            }
        }
        return(c(lower,upper))
    }
    colnames(history) <- c("Step", "Lower", "Upper")
    delta <- step.size/(1+nrow(history)/100)
    out$interval <- extract.interval(history, status, delta=delta, roll.width=roll.width)
    ##out$interval <- extract.interval(history, status, roll.width=roll.width)
    out$estimate <- mean(out$interval)
    out$history <- as_tibble(history)
    out$status <- status
    return(out)
}
