suppressMessages(library(tidyverse))
suppressMessages(library(bartMachine))
suppressMessages(library(glinternet))

t_stat_pval <- function(x.bar, s, n) {
    t_stat <- (x.bar) / (s / sqrt(n))
    2*pt(q=abs(t_stat), df=n-1, lower.tail=FALSE)
}

find_bart_interactions <- function(bart_machine_inter, idx_treat, idx_covar, num_replicates_for_avg = 25, delta_int=0.05, max_int=3) {
    p = ncol(bart_machine_inter$interaction_counts_avg)
    num.treat <- length(idx_treat)
    treatment = c()
    covariate = c()
    pvalue = c()
    for (i in idx_treat){
        for (j in idx_covar){
            treatment = c(treatment,i)
            covariate = c(covariate,j)
            x.bar = pmax(bart_machine_inter$interaction_counts_avg[i,j], bart_machine_inter$interaction_counts_avg[j,i])
            s = pmax(bart_machine_inter$interaction_counts_sd[i,j], bart_machine_inter$interaction_counts_sd[j,i])
            pvalue.new = t_stat_pval(x.bar, s, num_replicates_for_avg)
            pvalue = c(pvalue, pvalue.new)
        }
    }
    df = as_tibble(cbind(treatment,covariate,pvalue))
    df <- df %>% filter(pvalue<=delta_int) %>% group_by(treatment) %>% top_n(-max_int, pvalue)
    if(nrow(df)>0) {
        df = arrange(df, pvalue)
    } else {
        df = tibble(treatment=NA, covariate=NA, pvalue=NA) %>% head(0)
    }
    return(df)
}

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

cluster_covariates_lasso <- function(Z0, Z1, X.swap, Y, max_int=2) {
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
    X.full <- cbind(Z0, Z1, X.swap)
    p.full <- ncol(X.full)
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p.full

    interactionCandidates <- idx_treat
    interactionPairs <- lapply(idx_treat, function(i) {return(cbind(i,idx_covar))})
    interactionPairs <- do.call("rbind", interactionPairs)
    cv.fit = glinternet.cv(X.full, Y, rep(1,ncol(X.full)), nLambda=20, interactionPairs=interactionPairs, nFolds=2, screenLimit=max_int, verbose=TRUE)
    coef.cv = coef(cv.fit, s="lambda.1se")
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
            group_by(treatment) %>% top_n(max_int, interaction)
    }
    return(lasso.inter.tot)
}


cluster_covariates_BART <- function(Z0, Z1, X.swap, Y, delta_int=0.05) {
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
    X.full <- cbind(Z0, Z1, X.swap)
    p.full <- ncol(X.full)
    idx_covar <- (p0+1):(p0+p1)
    idx_treat <- (p0+p1+1):p.full

    bart_machine = bartMachine(data.frame(X.full), Y)
    num_replicates_for_avg = 10
    bart_machine_inter = interaction_investigator(bart_machine, num_trees_bottleneck=10, num_replicates_for_avg = num_replicates_for_avg, num_var_plot = 50, plot=TRUE)
    bart_inter = find_bart_interactions(bart_machine_inter, idx_treat, idx_covar, num_replicates_for_avg = num_replicates_for_avg, delta_int=delta_int)
    ## If the data are knockoff-augmented, rename the treatments
    bart_inter <- bart_inter %>% mutate(knockoff=treatment>p0+2*p1, treatment=treatment-(p0+p1), treatment=ifelse(knockoff,treatment-p1,treatment))
    bart_inter_tot <- bart_inter %>% group_by(treatment,covariate) %>%
        summarise(pvalue=min(pvalue)) %>%
        ungroup() %>%
        arrange(treatment,covariate)
    return(bart_inter_tot)
}
