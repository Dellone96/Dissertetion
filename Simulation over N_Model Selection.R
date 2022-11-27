library(magrittr)
library(invgamma)
library(mvtnorm)
library(xtable)


N = 50 # number of replicates

# Preliminary functions

comparison = function(true_gam, est_gam){
  
  out = c(sum(true_gam == 1 & est_gam == 1)/sum(true_gam), #TP Rate
          sum(true_gam == 0 & est_gam == 0)/sum(true_gam == 0), #TN Rate
          sum(true_gam == 0 & est_gam == 1)/sum(true_gam == 0), #FP Rate
          sum(true_gam == 1 & est_gam == 0)/sum(true_gam)) #FN Rate
  
  return(out)
  
}


m_y = function(y, X, a, b, beta0, V){
  
  p = ncol(X)
  n = length(y)
  
  beta0 = beta0[1:p]
  V     = V[1:p,1:p]
  
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


sim_over_samples = function(n, p, N){
  
  # Input:
  
  # n sample size
  # p number of covariates
  # N number of replicates
  
  # Storage:
  
  Gam_true = matrix(NA, nrow = N, ncol = p)
  
  Gam_bayes = matrix(NA, nrow = N, ncol = p)
  
  Gam_best_subset = matrix(NA, nrow = N, ncol = p)
  Gam_best_forward = matrix(NA, nrow = N, ncol = p)
  Gam_best_backward = matrix(NA, nrow = N, ncol = p)
  Gam_lasso = matrix(NA, nrow = N, ncol = p)
  Gam_FMA = matrix(NA, nrow = N, ncol = p)
  
  Y_all = list()
  X_all = list()
  
  Cfr_rates = array(NA, c(4, 6, N))
  
  #######################
  ## Replicates over N ##
  #######################
  
  for(i in 1:N){
    
    set.seed(i)
    
    # Generate true regression coefficients 
    
    beta_true = runif(p, -1, 1)
    
    gam  = c(1, sample(c(0,1), p - 1, replace = TRUE))
    
    ## Store:
    
    Gam_true[i,] = gam
    
    X = cbind(1, rmvnorm(n = n, mean = rep(0, p - 1), sigma = diag(1, p - 1)))
    
    # Generate respose variable from model Y = (beta*gam)'x + epsilon, with epsilon ~ N(0,1) 
    
    y = X%*%(beta_true*gam) + rnorm(n)
    
    ## Store:
    
    X_all[[i]] = X
    Y_all[[i]] = y
    
    
    # bayesian LR
    
    #### Marginal Distribution ####
    
    beta0  = rep(0, p)
    a = 0.01
    b = 0.01
    V = diag(1,p)
    
    final = c()
    
    # combine all the model 
    
    for (q in 1:(p-1)) {
      all_q = combn(1:(p-1),q)
      for (k in 1:ncol(all_q)) {
        t = m_y(y = y, X = X[,c(1,all_q[,k] + 1)], a = a, b = b, beta0 = beta0, V = V)  
        names(t) = paste(all_q[,k] + 1, collapse = ',')
        
        final = c(final,t)
      }
      
    }
    
    post = exp(final)/sum(exp(final))
    
    c = mean(final)
    
    post_c = exp(final - c)/(sum(exp(final - c)))
    
    #### Best bayesian model ####
    
    best_bayes = which.max(post_c)
    
    #### Gam bayesian comparison ####
    
    gam_bayes = rep(0,p)
    
    gam_bayes[as.integer(unlist(strsplit (names(best_bayes), ",")))]= 1
    
    gam_bayes[1] = 1
    
    Gam_bayes[i,] = gam_bayes
    
    # Frequentist Setting ----
    
    # Apply frequentist linear regression and see what is the model with the highest R squared index
    
    # Create a 0-1 vector of length p called gam_lm which has element j equal to 1 if Xj is included in the model, 0 otherwise
    # Compare gam_lm with gam to see if they match [**]
    
    #### Best Subsect ####
    
    library(ISLR)
    library(leaps)
    
    X = as.data.frame(X)
    
    out_subset  = regsubsets(y ~ ., X[,-1], nvmax = p)
    reg_summary = summary(out_subset)
    
    best_subset = which.max(reg_summary$adjr2)
    
    best_subset_pos = which(reg_summary$outmat[best_subset,] == "*")
    
    gam_best_subset = rep(0,p)
    gam_best_subset[1] = 1
    
    gam_best_subset[best_subset_pos + 1] = 1
    
    Gam_best_subset[i,] = gam_best_subset
    
    #### Forward Stepwise ####
    
    out_fwd = regsubsets(y ~ ., X, nvmax = p, method = "forward")
    reg_summary_fwd = summary(out_fwd)
    
    best_forward = which.max(reg_summary_fwd$adjr2)
    
    best_forward_pos = which(reg_summary$outmat[best_forward,] == "*")
    
    gam_best_forward = rep(0,p)
    gam_best_forward[1] = 1
    
    gam_best_forward[best_forward_pos + 1] = 1
    
    Gam_best_forward[i,] = gam_best_forward
    
    #### Backward ####
    
    out_bwd = regsubsets(y ~ ., X, nvmax = p, method = "backward")
    reg_summary_bwd = summary(out_bwd)
    
    best_bwd = which.max(reg_summary_bwd$adjr2)
    
    best_bwd_pos = which(reg_summary$outmat[best_bwd,] == "*")
    
    gam_best_bwd = rep(0,p)
    gam_best_bwd[1] = 1
    
    gam_best_bwd[best_bwd_pos + 1] = 1
    
    Gam_best_backward[i,] = gam_best_bwd
    
    #### Lasso ####
    
    library(glmnet)
    
    grid = 10^seq(-3, 9, length = 100)
    
    X = as.matrix(X)
    cv.out = cv.glmnet(X, y, lambda = grid, alpha = 1)
    
    bestlam = cv.out$lambda.min
    
    lasso_best = glmnet(X, y, alpha = 1, lambda = bestlam, standardize = TRUE)
    lasso_pos = lasso_best$beta@i
    
    gam_lasso = rep(0, p)
    gam_lasso[1] = 1
    
    gam_lasso[lasso_pos + 1] = 1
    
    Gam_lasso[i,] = gam_lasso
    
    
    #### For frequentist BMA ####
    
    # for each possible model fit the linear regression, compute the AIC (or BIC) instead of the marginal likelihood
    
    AICs = c()
    
    for (q in 1:(p-1)) {
      all_q = combn(1:(p-1),q)
      for (k in 1:ncol(all_q)) {
        
        #out_tmp = lm(y ~ X[,c(1,all_q[,k] + 1)] - 1)
        
        out_tmp = glm(y ~ X[,c(1,all_q[,k] + 1)] - 1)
        
        t = out_tmp$aic
        
        names(t) = paste(all_q[,k] + 1, collapse = ',')
        
        AICs = c(AICs,t)
      }
      
    }
    
    
    freq_weights = exp(-AICs/2)/sum(exp(-AICs/2))
    
    ## posso usarli per fare un BMA (su Y)
    
    best_FMA = which.min(AICs)
    
    gam_FMA = rep(0,p)
    
    gam_FMA[as.integer(unlist(strsplit (names(best_FMA), ",")))]= 1
    
    gam_FMA[1] = 1
    
    Gam_FMA[i,] = gam_FMA
    
    
    # Methods Comparison ----
    
    # Cfr tutti i metodi
    
    # Ogni metodo da un gamma stimato che cfr con il gamma vero in termini di rates
    
    
    bayesian = comparison(true_gam = gam, est_gam = gam_bayes)
    
    best_subset = comparison(true_gam = gam, est_gam = gam_best_subset)
    
    forward = comparison(true_gam = gam, est_gam = gam_best_forward)
    
    backward = comparison(true_gam = gam, est_gam = gam_best_bwd)
    
    lasso = comparison(true_gam = gam, est_gam = gam_lasso)
    
    FMA = comparison(true_gam = gam, est_gam = gam_FMA)
    
    table_comparison = data.frame(bayesian, best_subset, forward, backward, lasso, FMA )
    row.names(table_comparison) = c('True Positive', 'True Negative', 'False Positive', 'False Negative')
    
    
    Cfr_rates[,,i] = as.matrix(table_comparison)
    
    
  }
  
  return(out = list(Gam_true = Gam_true,
                    Gam_bayes = Gam_bayes,
                    Gam_best_subset = Gam_best_subset,
                    Gam_best_forward = Gam_best_forward,
                    Gam_best_backward = Gam_best_backward,
                    Gam_lasso = Gam_lasso,
                    Gam_FMA = Gam_FMA,
                    Y_all = Y_all,
                    X_all = X_all,
                    Cfr_rates = Cfr_rates))
  
}

### p = 5 ###

out_sim_n50_p5 = sim_over_samples(n = 50, p = 5, N = 50)
out_sim_n100_p5 = sim_over_samples(n = 100, p = 5, N = 50)
out_sim_n200_p5 = sim_over_samples(n = 200, p = 5, N = 50)
out_sim_n500_p5 = sim_over_samples(n = 500, p = 5, N = 50)
out_sim_n1000_p5 = sim_over_samples(n = 1000, p = 5, N = 50)

### p = 10 ###

out_sim_n50_p10 = sim_over_samples(n = 50, p = 10, N = 50)
out_sim_n100_p10 = sim_over_samples(n = 100, p = 10, N = 50)
out_sim_n200_p10 = sim_over_samples(n = 200, p = 10, N = 50)
out_sim_n500_p10 = sim_over_samples(n = 500, p = 10, N = 50)
out_sim_n1000_p10 = sim_over_samples(n = 1000, p = 10, N = 50)

### p = 15 ###

out_sim_n50_p15 = sim_over_samples(n = 50, p = 15, N = 50)
out_sim_n100_p15 = sim_over_samples(n = 100, p = 15, N = 50)
out_sim_n200_p15 = sim_over_samples(n = 200, p = 15, N = 50)
out_sim_n500_p15 = sim_over_samples(n = 500, p = 15, N = 50)
out_sim_n1000_p15 = sim_over_samples(n = 1000, p = 15, N = 50)


# Check su normalizzazione rates!
# p=5
average_rates_n50_p5 = apply(out_sim_n50_p5$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n100_p5 = apply(out_sim_n100_p5$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n200_p5 = apply(out_sim_n200_p5$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n500_p5 = apply(out_sim_n500_p5$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n1000_p5 = apply(out_sim_n1000_p5$Cfr_rates, FUN = mean, c(1,2), na.rm = T)

#p=10
average_rates_n50_p10 = apply(out_sim_n50_p10$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n100_p10 = apply(out_sim_n100_p10$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n200_p10 = apply(out_sim_n200_p10$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n500_p10 = apply(out_sim_n500_p10$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n1000_p10 = apply(out_sim_n1000_p10$Cfr_rates, FUN = mean, c(1,2), na.rm = T)

#p=15
average_rates_n50_p15 = apply(out_sim_n50_p15$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n100_p15 = apply(out_sim_n100_p15$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n200_p15 = apply(out_sim_n200_p15$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n500_p15 = apply(out_sim_n500_p15$Cfr_rates, FUN = mean, c(1,2), na.rm = T)
average_rates_n1000_p15 = apply(out_sim_n1000_p15$Cfr_rates, FUN = mean, c(1,2), na.rm = T)


### Creazione tabelle###

#p=5

table_p5 = rbind(round(average_rates_n50_p5, 2), 
                 round(average_rates_n100_p5, 2),
                 round(average_rates_n200_p5, 2), 
                 round(average_rates_n500_p5, 2),
                 round(average_rates_n1000_p5, 2))

row.names(table_p5)=c('TPR', 'TNR', 'FPR', 'FNR', 
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR')
colnames(table_p5)=c('Bayesian', ' Best_Subset', 'Forward','Backward', 'Lasso', 'FMA')

xtable(table_p5)


#p=10

table_p10 = rbind(round(average_rates_n50_p10, 2), 
                  round(average_rates_n100_p10, 2),
                  round(average_rates_n200_p10, 2), 
                  round(average_rates_n500_p10, 2),
                  round(average_rates_n1000_p10, 2))

row.names(table_p10)=c('TPR', 'TNR', 'FPR', 'FNR', 
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR')
colnames(table_p10)=c('Bayesian', ' Best_Subset', 'Forward','Backward', 'Lasso', 'FMA')
xtable(table_p10)


#p=15

table_p15 = rbind(round(average_rates_n50_p15, 2), 
                  round(average_rates_n100_p15, 2),
                  round(average_rates_n200_p15, 2), 
                  round(average_rates_n500_p15, 2),
                  round(average_rates_n1000_p15, 2))

row.names(table_p15)=c('TPR', 'TNR', 'FPR', 'FNR', 
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR',
                       'TPR', 'TNR', 'FPR', 'FNR')
colnames(table_p15)=c('Bayesian', ' Best_Subset', 'Forward','Backward', 'Lasso', 'FMA')
xtable(table_p15)

           

#save.image("out_all-10.Rdata")

load("out_all1610.Rdata")

