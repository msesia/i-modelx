suppressMessages(library(tidyverse))
suppressMessages(library(glmnet))

partition_covariates <- function(Z, vars.split) {
    bitsToInt <- function(x) {
        packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
    }
    n = nrow(Z)
    ## Number of partitions
    K = length(vars.split)
    out = rep(-1, n)
    for(i in 1:n) {
        out[i] = 1 + bitsToInt(Z[i,vars.split])
    }
    return(out)
}

aLKF_analysis <- function(Y, X, Xk, Z.noint, Z.int, ite=NULL, family = "binomial", fdr.nominal=0.1, fdr.offset=0,
                         naive=FALSE, split=FALSE, KFglobal=FALSE, cross.prior=TRUE, num.int=2, verbose=TRUE) {

    colnames(X) <- paste("X", 1:ncol(X), sep="_")
    colnames(Xk) <- paste("Xk", 1:ncol(Xk), sep="_")
    colnames(Z.noint) <- paste("Z", 1:ncol(Z.noint), sep="_")
    colnames(Z.int) <- paste("Z", ncol(Z.noint) + (1:ncol(Z.int)), sep="_")

    if(naive) cross.prior = FALSE
    if(split) cross.prior = FALSE
    if(KFglobal) cross.prior = TRUE

    random.swap=TRUE
    if(naive) random.swap = FALSE
    if(KFglobal) random.swap = FALSE

    n <- length(Y)
    N <- ncol(X)

    ## Split the data into two subsets (if applicable)
    ## Randomly swap variables and knockoffs (if applicable)
    if(split) {
        fold.1 <- sort(sample(n,round(n/2),replace=FALSE))
        fold.2 <- setdiff(1:n, fold.1)
        V.swap.1 <- matrix(rep(0,n*N), n)
        V.swap.2 <- matrix(rbinom(n*N,1,1/2), n)
    } else {
        if(naive) {
            V.swap <- matrix(rep(0,n*N), n)
        } else {
            V.swap <- matrix(rbinom(n*N,1,1/2), n)
        }
        fold.1 <- 1:n
        fold.2 <- 1:n
        V.swap.1 <- V.swap
        V.swap.2 <- V.swap
    }

    ## Partition the covariate space using the lasso with interactions
    if(KFglobal) {
        partition <- tibble(variable=NA, covariate=NA, covariate.name=NA, importance=NA) %>% head(0)
        env.list <- lapply(1:N, function(j) rep(1,n))
    } else {
        part.res <- aLKF_partition_intlasso(Y[fold.1], X[fold.1,], Xk[fold.1,], Z.noint[fold.1,], Z.int[fold.1,], V.swap.1[fold.1,],
                                           family=family, num.int=num.int, verbose=verbose)
        partition <- part.res$partition
        env.list <- part.res$env.list
    }

    ## Compute statistics for each variable in each group
    stats <- aLKF_compute_stats(Y[fold.2], X[fold.2,], Xk[fold.2,], Z.noint[fold.2,], Z.int[fold.2,], V.swap.2[fold.2,], partition,
                               family=family, random.swap=random.swap, KFglobal=KFglobal, cross.prior=cross.prior, scale.variables=TRUE, verbose=verbose)

    ## Make the full list of tested hypotheses and check which of them are truly null
    hypotheses <- aLKF_define_hypotheses(partition, env.list, N, ite=ite)

    ## Apply the knockoff filter
    disc.res <- aLKF_filter_stats(stats, fdr.nominal=fdr.nominal, fdr.offset=fdr.offset)
    discoveries <- disc.res$discoveries
    groups <- disc.res$groups

    ## Make the discoveries intellegible
    discoveries <- aLKF_parse_discoveries(discoveries, partition, N)

    ## Return results
    out <- c()
    out$discoveries <- discoveries
    out$hypotheses <- hypotheses
    return(out)
}

aLKF_define_hypotheses <- function(partition, groups, N, ite=NULL) {
    vars.split <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate})
    vars.split.names <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate.name})

    V.list <- lapply(1:N, function(j) {
        num_env = max(groups[[j]])
        w.j <- rep(j, num_env)
        return(w.j)
    })
    E.list <- lapply(1:N, function(j) {
        num_env = max(groups[[j]])
        w.j <- 1:num_env
        return(w.j)
    })

    hyp <- tibble(Variable=unlist(V.list))
    hyp$Splits <- sapply(1:nrow(hyp), function(i) {length(vars.split[[hyp$Variable[i]]])})
    hyp$Variables <- sapply(1:nrow(hyp), function(i) {
        paste(vars.split.names[[hyp$Variable[i]]], collapse=", ")
    })
    hyp$Partition=unlist(E.list)
    hyp$Label <- sapply(1:nrow(hyp), function(i) {
        if(length(vars.split[[hyp$Variable[i]]])>0) {
            idx.var <- sort(vars.split[[hyp$Variable[i]]])
            tmp.1 <- vars.split.names[[hyp$Variable[i]]]
            tmp.2 <- binaryLogic::as.binary(hyp$Partition[i]-1,n=hyp$Splits[i])
            tmp.3 <- paste(tmp.1, tmp.2,sep=" : ")
            out <- paste(tmp.3, collapse=";  ")
        } else {
            out <- "None"
        }
        return(out)
    })
    hyp <- hyp %>% select(-Splits)

     ## Check whether hypotheses are null
    hyp$Null <- NA
    hyp$Causal.prop <- NA
    hyp$Beta.sd <- NA
    hyp$n.sub <- NA
    for(i in 1:nrow(hyp)) {
        j <- hyp$Variable[i]
        idx <- which(groups[[j]]==hyp$Partition[i])
        if(is.null(ite)) {
            hyp$Null[i] <- NA
            hyp$Causal.prop[i] <- NA
            hyp$Beta.sd[i] <- NA
        } else {
            hyp$Null[i] = all(ite[idx,j]==0)
            hyp$Causal.prop[i] = mean(ite[idx,j]!=0)
            hyp$Causal.prop[i] <- ifelse(is.nan(hyp$Causal.prop[i]),NA,hyp$Causal.prop[i])
            hyp$Beta.sd[i] = sd(ite[idx,j])
        }
        hyp$n.sub[i] <- length(idx)
    }

    return(hyp)
}


aLKF_parse_discoveries <- function(discoveries, partition, N) {
    ## List of binary variable splits for each variable
    vars.split <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate})
    vars.split.names <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate.name})

    S.vec <- discoveries$Variable
    ## Make the results intellegible
    if(nrow(discoveries)>0) {
        ##discoveries <- tibble(Variable=V.vec[S.vec], Partition=E.vec[S.vec], W=W.vec[S.vec])
        discoveries$Splits <- sapply(1:nrow(discoveries), function(i) {length(vars.split[[discoveries$Variable[i]]])})
        discoveries$Label <- sapply(1:nrow(discoveries), function(i) {
            if(length(vars.split[[discoveries$Variable[i]]])>0) {
                idx.var <- vars.split[[discoveries$Variable[i]]]
                tmp.1 <- vars.split.names[[discoveries$Variable[i]]]
                tmp.2 <- binaryLogic::as.binary(discoveries$Partition[i]-1,n=discoveries$Splits[i])
                tmp.3 <- paste(tmp.1, tmp.2,sep=" : ")
                out <- paste(tmp.3, collapse=";  ")
            } else {
                out <- "None"
            }
            return(out)
        })
        discoveries <- discoveries %>% select(-Splits)
    } else {
        discoveries <- tibble(Offset=NA, Variable=NA, Partition=NA, Label=NA, W=NA) %>% head(0)
    }
    discoveries <- discoveries %>% arrange(Variable, desc(W))
    return(discoveries)
}

aLKF_partition_intlasso <- function(Y, X, Xk, Z.noint, Z.int, V.swap, num.int=2, family="gaussian", scale.variables=FALSE, verbose=TRUE) {

    if(verbose) {
        cat(sprintf("Looking for interactions with the lasso...\n"))
    }
    stopifnot(ncol(X)==ncol(Xk))
    N <- ncol(X)

    if(num.int==0) {
        cat(sprintf("Skipping interactions.\n"))
        out <- tibble(variable=1:N, covariate=NA, covariate.name=NA, importance=NA) %>% head(0)
        return(out)
    }

    ## Create data frame with all data
    data.tmp <- data.frame(cbind(Y,X,Xk,Z.noint, Z.int))
    colnames(data.tmp)[1] <- "Y"
    idx_treat <- 1+(1:(2*N))
    idx_covar.noint <- (max(idx_treat)+1):(max(idx_treat)+2)
    idx_covar.int <- (max(idx_covar.noint)+1):ncol(data.tmp)

    ## Make list of interactions
    interactionCandidates <- idx_treat
    interactionPairs <- lapply(idx_treat, function(i) {return(cbind(i,idx_covar.int))})
    interactionPairs <- as_tibble(do.call("rbind", interactionPairs))
    df.int <- interactionPairs %>% mutate(Variable = colnames(data.tmp)[i],
                                          Covariate = paste("`",colnames(data.tmp)[idx_covar.int],"`",sep="")) %>%
        select(-i, -idx_covar.int)

    ## Define the formula
    str.treat.1 <- paste(paste("X",1:N,sep="_"), collapse=" + ")
    str.treat.2 <- paste(paste("Xk",1:N,sep="_"), collapse=" + ")
    str.treat <- sprintf("%s + %s", str.treat.1, str.treat.2)
    str.noint <- paste(colnames(data.tmp)[idx_covar.noint], collapse=" + ")
    str.int <- paste(df.int$Variable, "*", df.int$Covariate, collapse=" + ")
    str.frmla <- sprintf("Y ~ %s + %s + %s", str.treat, str.noint, str.int)
    if(verbose) {
        cat(sprintf("Lasso formula:\n"))
        print(str.frmla)
    }

    ## Create the data matrix
    frmla <- as.formula(str.frmla)
    X.design <- model.matrix(frmla, data.tmp, keep.order = TRUE)[,-1]

    ## Create cloaked data set by swapping variables-knockoffs
    X.cloaked <- knockoff.random.swap(X.design[,1:N], X.design[,(N+1):(2*N)], V.swap)
    X.design.others <- X.design[,(2*N+1):ncol(X.design)]
    X.design.cloaked <- cbind(X.cloaked, X.design.others)

    ## Fit the lasso with interactions
    if(scale.variables) {
        cv.fit.int <- cv.glmnet(safe.scale(X.design.cloaked), Y, family=family, alpha=1)
    } else {
        cv.fit.int <- cv.glmnet(X.design.cloaked, Y, family=family, alpha=1)
    }
    ##plot(cv.fit.int)
    beta.hat.int <- coef(cv.fit.int, s="lambda.min")

    ## Parse the interaction model
    df <- tibble()
    for(variable in 1:N) {
        treat.label <- sprintf("X_%s:", variable)
        idx.int <- which(grepl(treat.label, rownames(beta.hat.int)))
        for(j in idx.int) {
            var.name <- strsplit(rownames(beta.hat.int)[j], ":", fixed=TRUE)[[1]][2]
            var.name <- str_replace_all(var.name, "`", "")
            covar.id <- which(colnames(Z.int)==var.name)
            coef.value <- beta.hat.int[j]
            df <- rbind(df, tibble(variable=variable, covariate=covar.id,
                                   covariate.name=var.name, importance=coef.value, knockoff=FALSE))
        }

        treat.label <- sprintf(":X_%s", variable)
        idx.int <- which(grepl(treat.label, rownames(beta.hat.int)))
        for(j in idx.int) {
            var.name <- strsplit(rownames(beta.hat.int)[j], ":", fixed=TRUE)[[1]][1]
            var.name <- str_replace_all(var.name, "`", "")
            covar.id <- which(colnames(Z.int)==var.name)
            coef.value <- beta.hat.int[j]
            df <- rbind(df, tibble(variable=variable, covariate=covar.id,
                                   covariate.name=var.name, importance=coef.value, knockoff=FALSE))
        }

        treat.label <- sprintf("Xk_%s:", variable)
        idx.int <- which(grepl(treat.label, rownames(beta.hat.int)))
        for(j in idx.int) {
            var.name <- strsplit(rownames(beta.hat.int)[j], ":", fixed=TRUE)[[1]][2]
            var.name <- str_replace_all(var.name, "`", "")
            covar.id <- which(colnames(Z.int)==var.name)
            coef.value <- beta.hat.int[j]
            df <- rbind(df, tibble(variable=variable, covariate=covar.id,
                                   covariate.name=var.name, importance=coef.value, knockoff=TRUE))
        }

        treat.label <- sprintf(":Xk_%s", variable)
        idx.int <- which(grepl(treat.label, rownames(beta.hat.int)))
        for(j in idx.int) {
            var.name <- strsplit(rownames(beta.hat.int)[j], ":", fixed=TRUE)[[1]][1]
            var.name <- str_replace_all(var.name, "`", "")
            covar.id <- which(colnames(Z.int)==var.name)
            coef.value <- beta.hat.int[j]
            df <- rbind(df, tibble(variable=variable, covariate=covar.id,
                                   covariate.name=var.name, importance=coef.value, knockoff=TRUE))
        }
    }
    inter.res <- df %>%
        group_by(variable, covariate, covariate.name) %>%
        summarise(importance=sum(abs(importance))) %>%
        ungroup() %>%
        arrange(variable, desc(importance))
    if(verbose) {
        cat(sprintf("List of interactions detected by the lasso:\n"))
        print(inter.res)
    }

    ## Extract the top interactions
    partition <- inter.res %>%
        filter(importance>0) %>%
        group_by(variable) %>%
        top_n(pmin(n(),num.int), importance)

    if(verbose) {
        cat(sprintf("List of selected interactions:\n"))
        print(inter.top)
    }

    ## List of binary variable splits for each variable
    vars.split <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate})
    ## List of individual clusters (environments) for each variable
    env.list <- lapply(1:N, function(j) { partition_covariates(Z.int, vars.split[[j]])})

    ## Return results
    out <- c()
    out$partition <- partition
    out$vars.split <- vars.split
    out$env.list <- env.list
    return(out)
}

aLKF_filter_stats <- function(stats, fdr.nominal=0.1, fdr.offset=NULL) {
    ## List of binary variable splits for each variable
    vars.split <- stats$vars.split
    ## List of individual clusters (environments) for each variable
    env.list <- stats$env.list

    ## Vectorize the statistics
    W.vec <- unlist(stats$W)
    V.vec <- unlist(stats$V)
    E.vec <- unlist(stats$E)
    W.sorted = W.vec[order(-abs(W.vec))]
    ## Apply the knockoff filter with different offsets
    if(is.null(fdr.offset)) {
        fdr.offset.values <- c(0,0.5,1)
    } else {
        fdr.offset.values <- c(fdr.offset)
    }
    discoveries <- tibble()
    for(fdr.offset in fdr.offset.values) {
        thres = knockoff.threshold(W.vec, fdr=fdr.nominal, offset=fdr.offset)
        S.vec <- which(W.vec >= thres)
        ##cat(sprintf("Made %d discoveries:\n", length(S.vec)))
        ##print(S.vec)

        ## Make the results intellegible
        if(length(S.vec)>0) {
            disc <- tibble(Offset=fdr.offset, Variable=V.vec[S.vec], Partition=E.vec[S.vec], W=W.vec[S.vec])
            disc$Splits <- sapply(1:nrow(disc), function(i) {length(vars.split[[disc$Variable[i]]])})
            disc$Label <- sapply(1:nrow(disc), function(i) {
                if(length(vars.split[[disc$Variable[i]]])>0) {
                    idx.var <- vars.split[[disc$Variable[i]]]
                    tmp.1 <- colnames(Z.int)[idx.var]
                    tmp.2 <- binaryLogic::as.binary(disc$Partition[i]-1,n=disc$Splits[i])
                    tmp.3 <- paste(tmp.1, tmp.2,sep=" : ")
                    out <- paste(tmp.3, collapse=";  ")
                } else {
                    out <- "None"
                }
                return(out)
            })
            disc <- disc %>% select(-Splits)
        } else {
            disc <- tibble(Offset=fdr.offset, Variable=NA, Partition=NA, W=NA, Label=NA) %>% head(0)
        }
        discoveries <- rbind(discoveries, disc)
    }
    out <- c()
    out$discoveries <- discoveries
    out$groups <- env.list
    return(out)
}

aLKF_compute_stats <- function(Y, X, Xk, Z.noint, Z.int, V.swap, partition, family = "binomial",
                              random.swap=TRUE, KFglobal=FALSE, cross.prior=TRUE, scale.variables=TRUE, verbose=TRUE) {

    N <- ncol(X)

    ## List of binary variable splits for each variable
    vars.split <- lapply(1:N, function(j) { filter(partition,variable==j)$covariate})
    ## List of individual clusters (environments) for each variable
    env.list <- lapply(1:N, function(j) { partition_covariates(Z.int, vars.split[[j]])})

    ## Compute the statistics group by group
    X.covar <- cbind(Z.noint, Z.int)
    if(scale.variables) {
        X <- safe.scale(X)
        Xk <- safe.scale(Xk)
        X.covar <- safe.scale(X.covar)
    }

    if(KFglobal) {
        ## Compute the test statistics simultaneously for all groups, with the usual lasso approach
        data.all <- cbind(X, Xk, X.covar)
        if(verbose) cat(sprintf("Computing all statistics simultaneously... "))
        cv.fit <- cv.glmnet(data.all, Y, family = family, alpha=1)
        if(verbose) cat(sprintf("Done.\n"))
        beta.hat = coef(cv.fit, s="lambda.min")[-1][1:(2*N)]
        stats <- c()
        stats$W <- abs(beta.hat[1:N]) - abs(beta.hat[(1+N):(2*N)])
        stats$V <- 1:N
        stats$E <- rep(1,N)
    } else {
        ## Compute the test statistics one variable at a time, group by group
        stats <- compute_stats_by_group(X.covar, X, Xk,
                                        Y, V.swap, env.list, family=family,
                                        calibrate.gamma=TRUE,
                                        cross.prior=cross.prior,
                                        random.swap=random.swap, verbose=verbose, dfmax=500)
    }
    stats$vars.split <- vars.split
    stats$env.list <- env.list
    return(stats)
}

compute_stats_by_group <- function(Z, X, Xk, Y, V.swap, env.list, family = "gaussian", calibrate.gamma=TRUE,
                                   cross.prior=FALSE, random.swap=TRUE, verbose=FALSE, dfmax=500) {
    p.z = ncol(Z)
    N = ncol(X)
    n = nrow(X)

    ## Create cloaked data set by swapping variables-knockoffs
    X.Xk.cloaked <- knockoff.random.swap(X, Xk, V.swap)
    X.cloaked <- X.Xk.cloaked[,1:N]
    Xk.cloaked <- X.Xk.cloaked[,(N+1):(2*N)]
    data.cloaked = cbind(Z, X.Xk.cloaked)
    ## Fit the lasso on the cloaked data
    beta.hat.prior <- tryCatch({
         cv.fit <- cv.glmnet(data.cloaked, Y, family = family, nfolds=10, alpha=1)
         beta.hat.prior = coef(cv.fit, s="lambda.min")[-1]
    }, error = function(e){
        beta.hat.prior = rep(1, ncol(data.cloaked))
    })
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
                ## Create cloaked data set by swapping variables-knockoffs, everywhere except in the current variable-environment
                idx.env = which(env.list[[j]]==k)
                V.swap.2 <- V.swap
                V.swap.2[idx.env,j] = 0
                X.cloaked.2 = X * (1-V.swap.2) + Xk * V.swap.2
                Xk.cloaked.2 = X * V.swap.2 + Xk * (1-V.swap.2)
            } else {
                X.cloaked = X
                Xk.cloaked = Xk
                X.cloaked.2 = X
                Xk.cloaked.2 = Xk
                idx.env = 1:n
            }
            if((cross.prior)) {
                if(length(idx.env)>0) {
                    ## Keep only the observations from this environment
                    data.X = cbind(Z, X.cloaked.2, Xk.cloaked.2)[idx.env,]
                    data.Y = Y[idx.env]
                    if(calibrate.gamma) {
                        ## Fit the lasso, tuning weight of prior
                        eval_gamma = function(gamma) {
                            penalty.z = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior.z)
                            penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
                            if((family=="binomial") && (min(table(data.Y))<2)) {
                                    return(c(0,1))
                            } else if((family=="gaussian") && (sd(data.Y)==0)) {
                                    return(c(0,1))
                            } else {
                                err <- tryCatch({
                                    cv_fit = cv.glmnet(data.X, data.Y, , family = family, nfolds=10, penalty.factor=c(penalty.z,rep(penalty,2)), alpha=1)
                                    idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
                                    c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
                                }, error = function(e){
                                    c(0,1)
                                })
                            }
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
                } else {
                    gamma = 0
                }
                ##if(verbose) cat(sprintf("Optimal prior weight: %.3f\n", gamma))

                ## Re-fit the lasso on environment k with optimally tuned prior
                penalty.z = 1*(1-gamma) + gamma * 1 / (0.05+beta.hat.prior.z)
                penalty = 1*(1-gamma) + gamma * 1 / (0.05+beta.hat.prior)
                if(length(idx.env)<2) {
                    w = 0
                } else {
                    if((family=="binomial") && (min(table(data.Y))<2)) {
                        w <- 0
                    } else if((family=="gaussian") && (sd(data.Y)==0)) {
                        w <- 0
                    } else {
                        w <- tryCatch({
                            cv.fit = cv.glmnet(data.X, data.Y, penalty.factor=c(penalty.z,rep(penalty,2)), alpha=1)
                            ## Extract the statistics for the variable of interest
                            beta.hat = coef(cv.fit, s="lambda.min")[-1][(p.z+1):(p.z+2*N)]
                            abs(beta.hat[j]) - abs(beta.hat[j+N])
                        }, error = function(e){
                            0
                        })
                    }
                }
            } else {
                if(min(table(Y))>=2) {
                    ## Fit the lasso on the cloaked data from all environments (without current variable)
                    data.cloaked = cbind(X.cloaked[,-j], Xk.cloaked[,-j],Z)
                    w <- tryCatch({
                        cv.fit.other <- cv.glmnet(data.cloaked, Y, family = family, nfolds=10, alpha=1)
                        beta.hat.other = coef(cv.fit.other, s="lambda.min")[-1]
                        ## Compute 'residuals' in current environment
                        offset <- data.cloaked[idx.env,] %*% beta.hat.other
                        ## Fit new model
                        lm.fit <- glm(Y[idx.env]~offset+X[idx.env,j]+Xk[idx.env,j],family=family)
                        beta.hat <- as.numeric(summary(lm.fit)$coefficients[-c(1,2),3])
                        w <- abs(beta.hat[1])-abs(beta.hat[2])
                        ifelse(is.na(w), 0, w)
                    }, error = function(e){
                        0
                    })
                } else {
                    w <- 0
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

knockoff.random.swap <- function(X, X.k, V.swap) {
    ## Create cloaked data set by swapping variables-knockoffs
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

evaluate_results <- function(discoveries, hypotheses) {
    fdr.offset.values <- unique(discoveries$Offset)
    ## Evaluate FDP
    discoveries <- discoveries %>% mutate(Offset=factor(Offset, fdr.offset.values))
    if(nrow(discoveries)>0) {
        discoveries.1 <- left_join(discoveries, hypotheses, by = c("Variable", "Partition", "Label"))
        res.tmp <- discoveries.1 %>% group_by(Offset) %>% summarise(FDP=mean(Null), Causal.prop=mean(Causal.prop,na.rm=T), Beta.sd=mean(Beta.sd,na.rm=T))
    } else {
        res.tmp <- tibble(Offset=fdr.offset.values[1], FDP=0, Causal.prop=NA, Beta.sd=NA) %>% mutate(Offset=factor(Offset, fdr.offset.values))
    }
    ## Evaluate power
    if(nrow(discoveries)>0) {
        discoveries.2 <- left_join(hypotheses, mutate(discoveries, Discovered=TRUE), by = c("Variable", "Partition", "Label")) %>%
            mutate(Discovered=ifelse(is.na(Discovered), FALSE, Discovered))
        res.tmp$Power <- discoveries.2 %>% filter(!Null) %>% summarise(Power=mean(Discovered)) %>% as.numeric()
        res.tmp$True <- discoveries.2 %>% filter(!Null) %>% summarise(Power=sum(Discovered)) %>% as.numeric()
        res.tmp$Discoveries <- discoveries.2 %>% summarise(Power=sum(Discovered)) %>% as.numeric()
    } else {
        res.tmp$Power <- 0
        res.tmp$True <- 0
        res.tmp$Discoveries <- 0
    }
    res.tmp <- res.tmp %>%
        complete(Offset) %>%
        mutate(FDP=ifelse(is.na(FDP), 0, FDP),
               True=ifelse(is.na(True), 0, True),
               Discoveries=ifelse(is.na(Discoveries), 0, Discoveries),
               Power=ifelse(is.na(Power), 0, Power))
    return(res.tmp)
}

safe.scale <- function(x) {
    apply(x, 2, function(y) {
        if(sd(y)>0) {
            return(scale(y))
        } else {
            return(y)
        }
    })
}

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
