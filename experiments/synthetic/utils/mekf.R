suppressMessages(library(glmnet))

knockoff.random.swap <- function(X, X.k, V.swap) {
    ## Create cloaked data set by swapping treatments-knockoffs
    n = nrow(X)
    X.aug <- cbind(X, X.k)
    X.out <- X.aug
    X.out.1 <- X*(1-V.swap)+X.k*(V.swap)
    colnames(X.out.1) <- colnames(X)
    X.out.2 <- X*(V.swap)+X.k*(1-V.swap)
    colnames(X.out.2) <- colnames(X.k)
    X.out <- cbind(X.out.1, X.out.2)
    return(X.out)
}

compute_stats_by_group <- function(Z, X, Xk, Y, V.swap, env.list, family = "gaussian", calibrate.gamma=TRUE,
                                   cross.prior=TRUE, residuals=FALSE, random.swap=TRUE, verbose=FALSE, dfmax=500, nfolds=10, alpha=1) {
    p.z = ncol(Z)
    N = ncol(X)
    n = nrow(X)
        ## ## Fit the lasso on all environments except k, then compute prior variable importance
        ## X_Xk_stack_mask = cbind(X_stack_mask, X_stack_k_mask)
        ## cv_fit = cv.glmnet(X_Xk_stack_mask[env_membership!=k,], y_stack[env_membership!=k])
        ## beta.hat.prior = coef(cv_fit, s="lambda.1se")[-1]
        ## beta.hat.prior = abs(beta.hat.prior[1:p])+abs(beta.hat.prior[(p+1):(2*p)])

        ## ## Fit the lasso on the data from environment k, tuning weight of prior
        ## X_env = cbind(Xs[[k]], X_ks[[k]])
        ## eval_gamma = function(gamma) {
        ##     penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
        ##     cv_fit = cv.glmnet(X_env, ys[[k]], penalty.factor=rep(penalty,2))
        ##     idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
        ##     err = c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
        ##     return(err)
        ## }
        ## gamma.seq = seq(0,1,length.out=10) # Sequence of prior weights to consider
        ## err.seq = sapply(gamma.seq, function(gamma) eval_gamma(gamma))

    ## Create cloaked data set by swapping treatments-knockoffs
    X.Xk.cloaked <- knockoff.random.swap(X, Xk, V.swap)
    X.cloaked <- X.Xk.cloaked[,1:N]
    Xk.cloaked <- X.Xk.cloaked[,(N+1):(2*N)]
    data.cloaked = cbind(Z, X.Xk.cloaked)
    ## Fit the lasso on the cloaked data
    cv.fit <- cv.glmnet(data.cloaked, Y, family = family, nfolds=nfolds, alpha=alpha)
    beta.hat.prior = coef(cv.fit, s="lambda.min")[-1]
    beta.hat.prior.z = abs(beta.hat.prior[1:p.z])
    beta.hat.prior = abs(beta.hat.prior[(p.z+1):(p.z+N)])+abs(beta.hat.prior[(p.z+N+1):(p.z+2*N)])

    ## Compute statistics for each variable and environment.
    W.list <- lapply(1:N, function(j) {
        num_env = max(env.list[[j]])
        w.j <- sapply(1:num_env, function(k) {
            if(verbose) {
                cat(sprintf("Computing statistics for variable %d in environment %d... ", j, k))
                flush.console()
            }
            if(random.swap) {
                ## Create cloaked data set by swapping treatments-knockoffs, everywhere except in the current variable-environment
                idx.env = which(env.list[[j]]==k)
                V.swap.2 <- V.swap
                V.swap.2[idx.env,j] = 0
                X.cloaked.2 = X * (1-V.swap.2) + Xk * V.swap.2
                Xk.cloaked.2 = X * V.swap.2 + Xk * (1-V.swap.2)
            } else {
                X.cloaked = X
                Xk.cloaked = Xk
            }
            if((cross.prior)) {
                ## Keep only the observations from this environment
                data.X = cbind(Z, X.cloaked.2, Xk.cloaked.2)[idx.env,]
                data.Y = Y[idx.env]
                if(calibrate.gamma) {
                    ## Fit the lasso, tuning weight of prior
                    eval_gamma = function(gamma) {
                        penalty.z = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior.z)
                        penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
                        cv_fit = cv.glmnet(data.X, data.Y, , family = family, nfolds=nfolds, penalty.factor=c(penalty.z,rep(penalty,2)), alpha=alpha)
                        idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
                        err = c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
                        return(err)
                    }
                    gamma.seq = seq(0,1,length.out=10) # Sequence of prior weights to consider
                    err.seq = sapply(gamma.seq, function(gamma) eval_gamma(gamma))
                    ## Find optimal prior weight
                    idx.best = which.min(colMeans(err.seq))
                    ## Make sure this is significantly better than gamma=0
                    if(err.seq[2,idx.best] >= mean(err.seq[,1])) {
                        idx.best  = 1
                    }
                    gamma = gamma.seq[idx.best]
                } else {
                    gamma = 0
                }
                ##if(verbose) cat(sprintf("Optimal prior weight: %.3f\n", gamma))

                ## Re-fit the lasso on environment k with optimally tuned prior
                penalty.z = 1*(1-gamma) + gamma * 1 / (0.05+beta.hat.prior.z)
                penalty = 1*(1-gamma) + gamma * 1 / (0.05+beta.hat.prior)
                if(length(data.Y)<2) {
                    w = 0
                } else {
                    if(sd(data.Y)==0) {
                        w = 0
                    } else {
                        cv.fit = cv.glmnet(data.X, data.Y, penalty.factor=c(penalty.z,rep(penalty,2)), alpha=alpha)
                        ## Extract the statistics for the treatment of interest
                        beta.hat = coef(cv.fit, s="lambda.min")[-1][(p.z+1):(p.z+2*N)]
                        w = abs(beta.hat[j]) - abs(beta.hat[j+N])
                    }
                }
            } else {
                if(residuals) {
                    if(length(idx.env)>0) {
                        ## Fit the lasso on the cloaked data from all environments (without current variable)
                        data.cloaked = cbind(X.cloaked[,-j], Xk.cloaked[,-j],Z)
                        cv.fit.other <- cv.glmnet(data.cloaked, Y, family = family, nfolds=nfolds, alpha=alpha)
                        beta.hat.other = coef(cv.fit.other, s="lambda.min")[-1]
                        ## Compute 'residuals' in current environment
                        offset <- data.cloaked[idx.env,] %*% beta.hat.other
                        ## Fit new model
                        lm.fit <- glm(Y[idx.env]~offset+X[idx.env,j]+Xk[idx.env,j],family=family)
                        beta.hat <- as.numeric(summary(lm.fit)$coefficients[-c(1,2),3])
                        w <- abs(beta.hat[1])-abs(beta.hat[2])
                        w <- ifelse(is.na(w), 0, w)
                    } else {
                        w <- 0
                    }
                } else {
                    data.X = cbind(Z, X.cloaked, Xk.cloaked)
                    data.Y = Y
                    cv.fit <- cv.glmnet(data.X, data.Y, alpha=alpha, family = family, nfolds=nfolds)
                    ## Extract the importance statistics for all treatments
                    beta.hat.all <- coef(cv.fit, s="lambda.min")[-1][(p.z+1):(p.z+2*N)]
                    ## Extract the statistics for the treatment of interest
                    w <- abs(beta.hat.all[j]) - abs(beta.hat.all[j+N])
                }
            }
            if(verbose) {
                cat(sprintf("%.3f\n", w))
                flush.console()
            }
            return(w)
        })
        return(w.j)
    })
    ## Return other info
    V.list <- lapply(1:N, function(j) {
        num_env = max(env.list[[j]])
        w.j <- rep(j, num_env)
        return(w.j)
    })
    E.list <- lapply(1:N, function(j) {
        num_env = max(env.list[[j]])
        w.j <- 1:num_env
        return(w.j)
    })
    out <- c()
    out$W <- W.list
    out$E <- E.list
    out$V <- V.list
    return(out)
}
