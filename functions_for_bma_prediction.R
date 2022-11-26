m_y = function(y, X, a, b, beta0, V){
  
  X = as.matrix(X)
  
  p = ncol(X)
  n = length(y)
  
  beta0 = beta0[1:p]
  V     = as.matrix(V[1:p,1:p])
  
  beta_tilde = solve(V+t(X)%*%X) %*% (t(X)%*%y + V%*%beta0)
  V_tilde = V + t(X)%*%X
  
  R = t(y) %*% y + t(beta0) %*% V %*% beta0 - t(beta_tilde) %*% (V+t(X) %*% X) %*% beta_tilde
  
  b_tilde = b + R
  
  a_tilde = a + n
  
  m = -n/2*log(2*pi) + 0.5*log(det(V)) - 0.5*log(det(V_tilde)) +
    0.5*a*log(0.5*b) - 0.5*a_tilde*log(0.5*b_tilde) +
    lgamma(0.5*a_tilde) - lgamma(0.5*a)
  
  return(as.integer(m))
  
}

## Ricorda: uso y_train e X_train per stimare il modello, ossia campionare dalla posterior dei parametri
## Uso X_test per il campionamento dalla predittiva

posterior_prediction_linear_conjugate = function(y_train, y_test, X_train, X_test, beta.0, V, a, b, S){
  
  X_train = as.matrix(X_train)
  X_test = as.matrix(X_test)
  
  ###########
  ## Input ##
  ###########
  
  # y     : response variable, (n,1) vector
  # X     : (n,p) design matrix, each column is one predictor, column one is the unit vector
  
  # beta.0, V : hyperparameters of the Multivariate Normal prior on beta
  # a, b      : hyperparameters of the Gamma prior on tau
  
  # S : number of MCMC iterations
  
  ############
  ## Output ##
  ############
  
  # out : a list with two elements:
  # beta_post   : (S,p) matrix, each row is a draw from the posterior of beta = (beta_1, ..., beta_p)
  # sigma2_post : (S,1) matrix, each row is a draw from the posterior of sigma2 = tau^{-1}
  
  
  library(mvtnorm)
  
  n = nrow(X_train)
  p = ncol(X_train)
  
  beta.0 = beta.0[1:p]
  V      = as.matrix(V[1:p,1:p])
  
  ## Compute posterior hyperparameters
  
  V.n    = V + t(X_train)%*%X_train
  beta.n = solve(V + t(X_train)%*%X_train)%*%(V%*%beta.0 + t(X_train)%*%y_train)
  
  a.n = a + n
  b.n = b + t(y_train)%*%y_train + beta.0%*%V%*%beta.0 - t(beta.n)%*%(V + t(X_train)%*%X_train)%*%beta.n
  
  ## To store samples drawn from the posterior of (beta, sigma2)
  
  beta_post   = matrix(NA, S, p)
  sigma2_post = matrix(NA, S, 1)
  
  Y_pred = matrix(NA, n, S)
  
  
  for(s in 1:S){
    
    #######################
    ## (1) Sample sigma2 ##
    #######################
    
    tau  = rgamma(1, a.n/2, b.n/2)
    sigma2 = 1/tau
    
    ##########################################
    ## (2) Sample beta conditionally on tau ##
    ##########################################
    
    beta = c(rmvnorm(1, beta.n, solve(V.n)/tau))
    
    #########################
    ## Store sampled draws ##
    #########################
    
    beta_post[s,]   = beta
    sigma2_post[s,] = sigma2
    
    
    ################
    ## Prediction ##
    ################
    
    y_hat = rnorm(n, X_test%*%beta, sigma2)
    
    Y_pred[,s] = y_hat
    
    
  }
  
  return(posterior = list(beta_post   = beta_post,
                          sigma2_post = sigma2_post,
                          Y_pred = Y_pred))
  
  
}
