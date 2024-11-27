## Computes empirical cross-prior statistics 
## Inputs:
##    Zs: a list of matrix Z from different environments
##    Xs: a list of matrix X from different environments
##    X_ks:  a list of knockoffs matrix X_k from different environments
##    ys:  a list of vector y from different environments
## Outputs:
##    A matrix of empirical cross-prior statistics
compute_stats_with_prior <- function(Zs, Xs, X_ks, Ys, family = "gaussian", verbose=FALSE, dfmax=500) {
    p = dim(Xs[[1]])[2]
    p.z = dim(Zs[[1]])[2]
    num_env = length(Xs)
    W_new_matrix = matrix(rep(0,num_env*p), ncol = num_env)

    for (k in 1:num_env){
        if(verbose) cat(sprintf("Computing statistics for environment %d...\n", k))

        X_stack_cloaked = NULL
        X_stack_k_cloaked = NULL
        y_stack = NULL
        env_membership = NULL
        for (k2 in 1:num_env) {
            ## Swap original-knockoffs in all environments (data cloaking)
            n = dim(Xs[[k2]])[1]
            M.swap = matrix(rbinom(n*p,1,1/2), n)
            X_cloaked = Xs[[k2]] * (1-M.swap) + X_ks[[k2]] * M.swap
            X_k_cloaked = Xs[[k2]] * M.swap + X_ks[[k2]] * (1-M.swap)
            X_stack_cloaked = rbind(X_stack_cloaked, X_cloaked)
            X_stack_k_cloaked = rbind(X_stack_k_cloaked, X_k_cloaked)
            y_stack = c(y_stack, Ys[[k2]])
            env_membership = c(env_membership, rep(k2, dim(Xs[[k2]])[1]))
        }

        ## Fit the lasso on all environments
        X_Xk_stack_cloaked = cbind(X_stack_cloaked, X_stack_k_cloaked)
        ##cv_fit = cv.glmnet(X_Xk_stack_cloaked[env_membership!=k,], y_stack[env_membership!=k], alpha=0)
        cv_fit = cv.glmnet(X_Xk_stack_cloaked, y_stack, alpha=0, family = family)
        beta.hat.prior = coef(cv_fit, s="lambda.min")[-1]
        beta.hat.prior = abs(beta.hat.prior[1:p])+abs(beta.hat.prior[(p+1):(2*p)])

        ## Fit the lasso on the data from environment k, tuning weight of prior
        X_env = cbind(Xs[[k]], X_ks[[k]])
        eval_gamma = function(gamma) {
            penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
            cv_fit = cv.glmnet(X_env, Ys[[k]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)
            idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
            err = c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
            return(err)
        }        
        gamma.seq = seq(0,1,length.out=10) # Sequence of prior weights to consider
        err.seq = sapply(gamma.seq, function(gamma) eval_gamma(gamma))

        if(FALSE) {
            tibble(gamma=gamma.seq, low=err.seq[1,], up=err.seq[2,]) %>%
                mutate(mid = (low+up)/2) %>%
                gather(low, up, mid, key="Limit", value="Err") %>%            
                ggplot(aes(x=gamma, y=Err, color=Limit)) +
                geom_point() +
                geom_line()
        }        
        
        ## Find optimal prior weight
        idx.best = which.min(colMeans(err.seq))
        ## Make sure this is significantly better than gamma=0
        if(err.seq[2,idx.best] >= mean(err.seq[,1])) {
            idx.best  = 1
        }
        gamma = gamma.seq[idx.best]
        if(verbose) cat(sprintf("Optimal prior weight: %.3f\n", gamma))
        
        ## Re-fit the lasso on environment k with optimally tuned prior
        penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
        cv_fit = cv.glmnet(X_env, Ys[[k]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)

        ## Extract importance measures
        beta.hat = coef(cv_fit, s="lambda.min")[-1][1:(2*p)]
        W1 = abs(beta.hat[1:p]) - abs(beta.hat[(p+1):(2*p)])
        W_new_matrix[,k] = W1
    }
    return(W_new_matrix)
}

################################################################
## These functions run the multi-environment knockoff analysis
################################################################
library(knockoff)

## Report discoveries from testing consistent conditional associations
## Inputs: 
##    W_matrix: a matrix of multi-environment knockoff statistics
##    q: FDR threshold
## Outputs:
##    Variables selected by the multi-environment knockoff filter
invariant_model = function(W_matrix, q = 0.1, offset=1){
    W_sign = 2*((apply(W_matrix, 1, min) > 0) - 0.5)
    W_prod = abs(apply(W_matrix, 1, prod))
    W_eff = W_sign * W_prod
    thres = knockoff.threshold(W_eff, fdr=q, offset=offset)
    selected = which(W_eff >= thres)
    return(selected)
}

## Find the sign of entries in a vector
## Inputs: 
##    vec: a vector
## Outputs:
##    The sign of the each entry: -1 for negative, +1 for positive, and +1/-1
##              with probability 0.5 for 0. 
sign_fun = function(vec){
    l = length(vec)
    return((rbinom(l,1,0.5)*(vec == 0) + (vec > 0))*2 - 1)
}

## Compute partial conjunction multi-environment p-values
## Inputs: 
##    vec: a vector of knockoff statistics
##    r: the number of nonnull environments required
##    randomness: randomness = 0 computes the p-values as in (16) in the paper without the Uj term
##                randomness = 1 computes the p-values as in (19) in the paper without the Uj term
##                randomness = 2 computes the p-values as in (16) in the paper with the Uj term
##                randomness = 3 computes the p-values as in (19) in the paper with the Uj term      
## Outputs:
##    Partial conjunction multi-environment p-values
p_value = function(vec, r, randomness){
    l = length(vec)
    minus = sum(vec < 0)
    num_zero = sum(vec == 0)
    
    if(randomness == 3){
        minus = minus + rbinom(1,num_zero,0.5)
        num_zero = 0
    }
    
    a0 = l - r + 1
    a = l - r + 1 - num_zero
    if (a > 0){
        discrete_p = pbinom(minus,a,0.5)
        discrete_p_l = pbinom(minus - 1,a,0.5)
    } else {
        discrete_p = 1
        discrete_p_l = 0
    }
    p_cont = runif(1,discrete_p_l,discrete_p)
    steps = pbinom(0:a0,a0,0.5)
    p_disc = steps[min(which(steps >= p_cont))]
    if(randomness == 0){
        return (discrete_p)
    }else if(randomness == 1){
        return (p_disc)
    }else if(randomness >= 2){
        return (p_cont)
    }
}

## Combine multi-environment knockoff statistics
## Inputs: 
##    vec: a vector of knockoff statistics
##    r: the number of nonnull environments required
## Outputs:
##    invariant statistics
combine_mag = function(vec, r){
    vec_abs = abs(vec)
    vec_sort = -sort(-vec_abs, partial = r)[1:r]
    return(prod(vec_sort))
}


## Report discoveries from testing partically consistent conditional associations
## Inputs: 
##    W_matrix: a matrix of multi-environment knockoff statistics
##    r: the number of nonnull environments required
##    q: FDR threshold
##    method: variable section method; "seqstep" stands for Selective Seqstep+, "accumulation" stands for Accumulation test
##    c: threshold c in Selective Seqstep+
##    randomness: randomness = 0 computes the p-values as in (16) in the paper without the Uj term
##                randomness = 1 computes the p-values as in (19) in the paper without the Uj term
##                randomness = 2 computes the p-values as in (16) in the paper with the Uj term
##                randomness = 3 computes the p-values as in (19) in the paper with the Uj term      
## Outputs:
##    Variables selected by the multi-environment knockoff filter
partial_conjunction = function(W_matrix, r, q = 0.1, method = "seqstep", c = 0.6, randomness = 2, offset=0){
    pvals = apply(W_matrix, 1, p_value, r = r, randomness = randomness)
    W_sign = 2*(pvals <= c) - 1
    W_mag = apply(W_matrix, 1, combine_mag, r = r)
    if (method == "seqstep"){
        q_tilde = q*(1-c)/c
        W_eff = W_sign * W_mag
        thres = knockoff.threshold(W_eff, fdr=q_tilde, offset=offset)
        selected = which(W_eff >= thres)
        return (selected)
    } else if(method == "accumulation"){
        hfun = create_HingeExp_function(C=2)
        pvals_sorted = pvals[order(-W_mag)]
        num_selecteda = AccumulationTest(pvals_sorted , hfun, alpha=q)
        selecteda = NULL
        if (num_selecteda > 0){
            selecteda = order(-W_mag)[1:num_selecteda]
        }
        return(sort(selecteda))
    }
    return (NULL)
}

##############################################################################
## These functions implement 
##   the Accumulation Test methods from the paper:
## Ang Li & Rina Foygel Barber,
##   "Accumulation tests for FDR control
##      in ordered hypothesis testing"
## Available from http://arxiv.org/abs/1505.07352
## (Several methods from other papers also implemented,
##        as noted below - see citations in paper)
##############################################################################

##############################################################################
## HingeExp method,
##    i.e. an accumulation test with the HingeExp function:
##       h(p) = C * log(1/(C*(1-p))) * (p>1-1/C)
##############################################################################
create_HingeExp_function = function(C=2){
	function(x){C*log(1/(C*(1-x)))*(x>1-1/C)}
}
HingeExp = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_HingeExp_function(C),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## ForwardStop method (G'Sell et al 2013),
##    i.e. an accumulation test with the ForwardStop function:
##       h(p) = log(1/(1-p))
##############################################################################
create_ForwardStop_function = function(){
	function(x){log(1/(1-x))}
}
ForwardStop = function(pvals,alpha=0.2,output_type='khat'){
	AccumulationTest(pvals,create_ForwardStop_function(),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## SeqStep method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##############################################################################
create_SeqStep_function = function(C=2){
	function(x){C*(x>1-1/C)}
}
SeqStep = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_SeqStep_function(C),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
####################################################################################


##############################################################################
## SeqStep+ method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##         & with the conservative correction
##              for estimating FDR
##############################################################################
SeqStepPlus = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_SeqStep_function(C),alpha=alpha,numerator_plus=C,denominator_plus=1,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## Accumulation test for a generic function "hfun"
##############################################################################
AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat',check_integrate_to_one=TRUE){
	
	# check for valid arguments
	check_inputs = CheckInputs_AccumulationTest(pvals,hfun,alpha,numerator_plus,denominator_plus,output_type,check_integrate_to_one)
	if(length(check_inputs)>0){
		stop(check_inputs)
	}
	
	# perform the test
	n=length(pvals)
	FDPest=(numerator_plus+cumsum(unlist(lapply(pvals,hfun))))/(denominator_plus+1:n)
	FDPest_vs_alpha=(FDPest%*%t(rep(1,length(alpha)))<=rep(1,n)%*%t(alpha))
	findlast=function(x){max(c(0,which(x)))}
	khat=apply(FDPest_vs_alpha,2,findlast)
	if(output_type=='khat'){		
		return(khat)
	}else{if(output_type=='FDPest'){
		return(FDPest)
	}else{
		output=list()
		output$FDPest=FDPest
		output$khat=khat
		return(output)
	}}
}
##############################################################################



##############################################################################
## Check inputs for AccumulationTest
##############################################################################
CheckInputs_AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat', check_integrate_to_one){
	# check_integrate_to_one should be logical
	if(!is.logical(check_integrate_to_one)){
		return('check_integrate_to_one must be logical')
	}
	# check that pvals and alpha are each sequences with values in [0,1]
	if(!is.numeric(pvals) || !is.vector(pvals) || min(pvals)<0 || max(pvals)>1){
		return('pvals must be a number or numeric vector with values in [0,1]')
	}
	n=length(pvals)
	
	if(!is.numeric(alpha) || !is.vector(alpha) || min(alpha)<0 || max(alpha)>1){
		return('alpha must be a number or numeric vector with values in [0,1]')
	}	
	
	# check that hfun is a function that gives a nonnegative value for each pvalue
	if(!is.function(hfun)){
		return('hfun must be a function')
	}
	if(!is.numeric(try(hfun(pvals),silent=TRUE)) || any(is.na(try(hfun(pvals),silent=TRUE)))){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(length(hfun(pvals))!=length(pvals)){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(any(hfun(pvals)<0)){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(check_integrate_to_one){
		if(abs(integrate(hfun,0,1)$value-1)>1e-2){
			return('The function hfun must have expectation 1 when applied to a uniform variable p~Uniform[0,1] (set check_integrate_to_one=FALSE to override)')
		}
	}
	
	# check that numerator_plus and denominator_plus are numbers
	if(!is.numeric(numerator_plus) || !is.numeric(denominator_plus) || length(numerator_plus)!=1 || length(denominator_plus)!=1){
		return('numerator_plus and denominator_plus must each be a scalar')
	}
	
	# check that output_type is in {'khat', 'FDPest', 'both'}
	if(!is.element(output_type,c('khat','FDPest','both'))){
		return('Invalid output type: choose from "khat", "FDPest", or "both"')
	}
	
	return(NULL)
}
##############################################################################
