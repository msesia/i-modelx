suppressMessages(library(tidyverse))
suppressMessages(library(grf))
suppressMessages(library(glmnet))
suppressMessages(library(glinternet))

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

cluster_covariates_lasso <- function(Z0, Z1, X.swap, Y, fold.1, fold.2, num_cluster=4) {
    #################################
    ## Z0: passive covariates
    ## Z1: interaction covariates
    ## X.swap: treatment variables (possibly augmented with knockoffs)
    ## Y: outcome
    #################################
    ## Combine data
    p0 <- ncol(Z0)
    p1 <- ncol(Z1)
    p.x <- ncol(X.swap)
    p.treat <- p.x/2
    X.full <- cbind(Z0, Z1, X.swap)
    p.full <- ncol(X.full)
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p.full

    interactionCandidates <- idx_treat
    interactionPairs <- lapply(idx_treat, function(i) {return(cbind(i,idx_covar))})
    interactionPairs <- do.call("rbind", interactionPairs)
    screenLimit=round(log(num_cluster)/log(2))
    cv.fit = glinternet.cv(X.full[fold.1,], Y[fold.1], rep(1,ncol(X.full)), nLambda=20, interactionPairs=interactionPairs, nFolds=2, screenLimit=screenLimit, verbose=TRUE)
    coef.cv = coef(cv.fit, s="lambda.min")
    if(is.null(coef.cv$interactions$contcont)) {
        lasso.inter.tot = tibble(treatment=NA, covariate=NA, interaction=NA) %>% head(0)
    } else {
        lasso.inter = tibble(treatment=coef.cv$interactions$contcont[,2], covariate=coef.cv$interactions$contcont[,1])
        lasso.inter$interaction = unlist(coef.cv$interactionsCoef$contcont)
        lasso.inter = lasso.inter %>% arrange(treatment, covariate)
        ## If the data are knockoff-augmented, rename the treatments
        lasso.inter <- lasso.inter %>% mutate(knockoff=treatment>p0+2*p1, treatment=treatment-(p0+p1), treatment=ifelse(knockoff,treatment-p1,treatment))
        lasso.inter.tot <- lasso.inter %>% group_by(treatment,covariate) %>%
            summarise(interaction=sum(abs(interaction))) %>%
            ungroup() %>%
            arrange(treatment,covariate) %>%
            group_by(treatment) %>% top_n(screenLimit, interaction)            
    }
    ## List of binary variable splits for each treatment
    vars.split <- lapply(1:p1, function(j) {
        filter(lasso.inter.tot,treatment==j)$covariate
    })
    ## List of individual clusters (environments) for each variable
    Z <- cbind(Z0, Z1)
    env.list <- lapply(1:p1, function(j) { partition_covariates(Z[fold.2,], vars.split[[j]])})
        
    if(FALSE) {

        clusters = env.list[[2]]
        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][1] ])})
        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][2] ])})
        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][1] ]*Z1[which(clusters==k),interaction_design[[j]][2] ])})

        Z1[which(beta_e[idx_treat[j],]!=0), interaction_design[[j]][1] ]
        Z1[which(beta_e[idx_treat[j],]!=0), interaction_design[[j]][2] ]

    }
    return(env.list)
}    

cluster_covariates_ite <- function(j, Z0, Z1, X.swap, Y, fold.1, fold.2, num_cluster=4) {
    #################################
    ## Z0: passive covariates
    ## Z1: interaction covariates
    ## X.swap: treatment variables (possibly augmented with knockoffs)
    ## Y: outcome
    #################################
    ## Combine data
    p0 <- ncol(Z0)
    p1 <- ncol(Z1)
    p.x <- ncol(X.swap)
    p.treat <- p.x/2
    X.full <- cbind(Z0, Z1, X.swap)
    p.full <- ncol(X.full)
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p.full

    ## Fit causal random forest
    estimate_ite <- function(var.idx) {
        data.covar <- cbind(Z0,Z1)
        W <- X.swap[,var.idx]
        W.hat <- Z0[fold.1,1:p1][,j]
        tau.forest <- causal_forest(data.covar[fold.1,], Y[fold.1], W[fold.1], W.hat=W.hat, num.trees=100, compute.oob.predictions=FALSE)
        ## Estimate ITE (for second fold)
        ite.fold.1 <- predict(tau.forest, data.covar[fold.1,])[,1]
        ite.fold.2 <- predict(tau.forest, data.covar[fold.2,])[,1]
        ite <- c()
        ite$fold.1 <- ite.fold.1
        ite$fold.2 <- ite.fold.2
        return(ite)
    }
    ite.1 <- estimate_ite(j)
    ite.2 <- estimate_ite(j+p.treat)
    ite.fold.1 <- ite.1$fold.1 + ite.2$fold.1
    ite.fold.2 <- ite.1$fold.2 + ite.2$fold.2
    ## Cluster the ITE
    quantile.levels <- seq(1/num_cluster, (num_cluster-1)/num_cluster, length.out=num_cluster-1)
    quantiles.fold.1 <- quantile(ite.fold.1, quantile.levels)
    cut.points <- c(-Inf, quantiles.fold.1, Inf)
    clusters <- as.numeric(cut(ite.fold.2, cut.points, labels=1:num_cluster)   )
    
    if(FALSE) {

        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][1] ])})
        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][2] ])})
        sapply(1:4, function(k) {mean(Z1[which(clusters==k),interaction_design[[j]][1] ]*Z1[which(clusters==k),interaction_design[[j]][2] ])})

        Z1[which(beta_e[idx_treat[j],]!=0), interaction_design[[j]][1] ]
        Z1[which(beta_e[idx_treat[j],]!=0), interaction_design[[j]][2] ]

    }
    return(clusters)
}
