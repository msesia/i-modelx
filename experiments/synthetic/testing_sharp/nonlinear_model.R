generate_tree = function(n, p, p1) {
  ## Problem parameters
  signal_mean = 10  # average effect size for causal treatments
  signal_std = 0   # standard deviation of effect size for causal treatments
  batch <- 1
  batch.size <- 1
  
  seed = 2022
  seed.model = 2022
  
  ## Parameters for P(X)
  rho = 0.5
  p0 = p - p1           # Number of covariates
  
  ## Parameters for P(Y|X)
  p_causal = round(p1/2)        # Number of causal treatments
  p_causal_covar = round(p0/2)  # Number of causal covariates
  sample_X = function(n, p, p1) {
    ## Number of non-treatment variables
    p0 = p - p1
    stopifnot(p1 <= p0)
    ## Generate half of the variables from a multivariate normal distribution
    mu = rep(0,p0)
    Sigma = toeplitz(rho^(0:(p0-1)))
    X1 = matrix(rnorm(n*p0),n) %*% chol(Sigma)
    X1 = pnorm(X1-1)
    # Generate the second half of the variables from Bernoulli conditional on the first half
    U2 = matrix(runif(p1*n), n, p1)
    X2 = (U2 <= X1[,1:p1])
    # Combine the variables
    X = cbind(X2, X1)
    return(X)
  }
  
  sample_effects = function(p, p1, p_causal, p_causal_covar, signal_mean, signal_std) {
    p0 = p - p1
    nonzero_treat = sort(sample(p1, p_causal))
    nonzero_covar = sort(p1+sample(p0, p_causal_covar))
    nonzero = sort(c(nonzero_treat, nonzero_covar))
    # Flip signs for covariates
    signs = c(rep(1,p1), 2*rbinom(p0,1,0.5)-1)
    beta_center = signal_mean * (1:p %in% nonzero) * signs
    beta_center = matrix(rep(beta_center,each=n), n, p)
    beta_matrix = beta_center + matrix(rnorm(n*p, mean=0, sd=signal_std), n, p)
    beta_matrix[,-nonzero] = 0
    beta_matrix = t(beta_matrix)
    return(beta_matrix)
  }
  
  sample_Y_linear = function(X, beta, interaction_strength=0) {
    n = nrow(X)
    p = ncol(X)
    y = rep(0,n)
    for(i in 1:n) {
      y[i] = X[i,,drop=FALSE] %*% beta[,i,drop=FALSE]
    }
    if(interaction_strength>0) {
      beta_avg = rowMeans(beta)[(p1+1):p]
      beta_support = sort(which(abs(beta_avg) > 0))
      pairs = lapply(seq(1,5), function(tmp) p1 + sample(beta_support, 2, replace=FALSE))
      for(pair in pairs) {
        y = y + interaction_strength * (X[,pair[1]] * X[,pair[2]] - mean( X[,pair[1]] * X[,pair[2]]) )
      }
    }
    y = y + rnorm(n)
    return(y)
  }
  
  ## Set random seed for the model
  set.seed(seed.model)
  
  ## Sample the effect matrix
  beta = sample_effects(p, p1, p_causal, p_causal_covar, signal_mean, signal_std)
  rowMeans(beta)
  
  ## Set random seed for the experiment
  set.seed(seed)
  
  ## Sample the data
  X = sample_X(n, p, p1)
  Y = sample_Y_linear(X, beta, interaction_strength=5)
  
  library(rpart)
  library(rpart.plot)
  
  df = as_tibble(X[,seq(p1+1,p)]) %>% mutate(Y=Y)
  tree = rpart(Y~., df)
  
  #identify best cp value to use
  best <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  
  #produce a pruned tree based on the best cp value
  pruned_tree <- prune(tree, cp=best)
  
  if(TRUE) {
    prp(tree,
        faclen=0, #use full names for factor labels
        extra=1, #display number of obs. for each terminal node
        roundint=F, #don't round to integers in output
        digits=5) #display 5 decimal places in output
    #plot the pruned tree
    prp(pruned_tree,
        faclen=0, #use full names for factor labels
        extra=1, #display number of obs. for each terminal node
        roundint=F, #don't round to integers in output
        digits=5) #display 5 decimal places in output
  }  
  return(tree)
}
