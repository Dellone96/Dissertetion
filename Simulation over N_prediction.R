library(magrittr)
library(invgamma)
library(mvtnorm)
library(xtable)
library(ISLR)
library(leaps)
library(glmnet)

#N = 5 # number of replicates

#n = 100
#p = 5

# Preliminary functions

source("functions_for_bma_prediction.R")

sim_over_samples_pred = function(n, p, N){
  
  # Input:
  
  # n sample size
  # p number of covariates
  # N number of replicates
  
  # Storage:
  
  Gam_true = matrix(NA, nrow = N, ncol = p)
  
  MSE_best_subset = matrix(NA, nrow = N, ncol = 1)
  MSE_forward = matrix(NA, nrow = N, ncol = 1)
  MSE_bwd = matrix(NA, nrow = N, ncol = 1)
  MSE_lasso = matrix(NA, nrow = N, ncol = 1)
  MSE_FMA = matrix(NA, nrow = N, ncol = 1)
  
  MSE_BMA = matrix(NA, nrow = N, ncol = 1)
  
  Y_all = list()
  X_all = list()
  
  #######################
  ## Replicates over N ##
  #######################
  
  for(i in 1:N){
    
    set.seed(i)
    
    #Generate coeff
    
    beta_true = runif(p, -1, 1)
    gam  = c(1, sample(c(0,1), p - 1, replace = TRUE))
    Gam_true[i,] = gam
    
    X = cbind(1, rmvnorm(n = n, mean = rep(0, p - 1), sigma = diag(1, p - 1)))
    
    # Generate response variable from model Y = (beta*gam)'x + epsilon, with epsilon ~ N(0,1) 
    
    y = X%*%(beta_true*gam) + rnorm(n)
    
    train = sample(1:n, n/2) # these observations will be used to fit the model
    test = setdiff(1:n, train) # difference between sets 1:n and train
    
    X_train = X[train,]
    X_test  = X[test,]
    
    y_train = y[train]
    y_test = y[test]
    
    ## Store:
    
    X_all[[i]] = X
    Y_all[[i]] = y
    
    
    ###Bayesian##
    
    beta.0 = beta0 = rep(0, p)
    a = 0.01
    b = 0.01
    V = diag(1,p)
    
    S = 500
    
    ## Repeat for all possible models corresponding to subsets of columns of X
    ## that is X[,c(1,all_q[,k] + 1)] inside the loop (as in the FMA with glm)
    
    # combine all the model 
    
    Y_pred_bma = c()
    
    final = c()
    
    ## Model with no covariates (just the intercept)
    
    t = m_y(y = y_train, X = X_train[,c(1), drop = FALSE], a = a, b = b, beta0 = beta0, V = V)  # again, m_y is for model fitting (then use train)
    names(t) = paste("0")
    
    out_k = posterior_prediction_linear_conjugate(y_train, y_test, X_train[,c(1), drop = FALSE], X_test[,c(1), drop = FALSE], beta.0, V, a, b, S)
    
    y_pred_k = rowMeans(out_k$Y_pred)
    
    Y_pred_bma = cbind(Y_pred_bma, y_pred_k)
    
    
    final = c(final,t)
    
    ##
    
    for (q in 1:(p-1)) {
      all_q = as.matrix(combn(1:(p-1),q))
      for (k in 1:ncol(all_q)) {
        t = m_y(y = y_train, X = X_train[,c(1,all_q[,k] + 1)], a = a, b = b, beta0 = beta0, V = V)  # again, m_y is for model fitting (then use train)
        names(t) = paste(all_q[,k] + 1, collapse = ',')
        
        out_k = posterior_prediction_linear_conjugate(y_train, y_test, X_train[,c(1,all_q[,k] + 1)], X_test[,c(1,all_q[,k] + 1)], beta.0, V, a, b, S)
        
        y_pred_k = rowMeans(out_k$Y_pred)
        
        Y_pred_bma = cbind(Y_pred_bma, y_pred_k)
        
        
        final = c(final,t)
      }
      
    }
    
    post = exp(final)/sum(exp(final))
    
    c = mean(final)
    
    post_c = exp(final - c)/(sum(exp(final - c)))
    
    y_hat_bma = Y_pred_bma%*%post_c
    
    mse_bma = mean((y_hat_bma - y_test)^2)
    MSE_BMA[i,] = mse_bma
    
    
    ###Best Subset###
    
    X = as.data.frame(X)
    
    out_subset  = regsubsets(y ~ ., X[,-1], nvmax = p)
    reg_summary = summary(out_subset)
    
    best_subset = which.max(reg_summary$adjr2)
    
    best_subset_pos = which(reg_summary$outmat[best_subset,] == "*")
    
    gam_best_subset = rep(0,p)
    gam_best_subset[1] = 1
    
    gam_best_subset[best_subset_pos + 1] = 1
    gam_best_subset
    
    ## Evaluate the accuracy by fitting the model on the train and then using the test
    
    
    out_best_subset = lm(y ~.-1, data= X[,which(gam_best_subset == 1)], subset = train)
    
    X_test = as.data.frame(X_test)
    
    y_hat_best_subset = predict(out_best_subset, newdata = X_test[,which(gam_best_subset == 1)])
    
    mse_best_subset = mean((y_test - y_hat_best_subset)^2)
    mse_best_subset
    
    
    MSE_best_subset[i,] = mse_best_subset
    
    ### Forward ###
    
    out_fwd = regsubsets(y ~ .-1, X, nvmax = p, method = "forward")
    reg_summary_fwd = summary(out_fwd)
    
    best_forward = which.max(reg_summary_fwd$adjr2)
    
    best_forward_pos = which(reg_summary$outmat[best_forward,] == "*")
    
    gam_best_forward = rep(0,p)
    gam_best_forward[1] = 1
    
    gam_best_forward[best_forward_pos + 1] = 1
    gam_best_forward
    
    out_best_forward = lm(y ~.-1, data= X[,which(gam_best_forward == 1)], subset = train)
    
    X_test = as.data.frame(X_test)
    
    y_hat_best_forward = predict(out_best_forward, newdata = X_test[,which(gam_best_forward == 1)])
    
    mse_best_forward = mean((y_test - y_hat_best_forward)^2)
    MSE_forward[i,] = mse_best_forward
    
    ### Backward ###
    
    out_bwd = regsubsets(y ~ .-1, X, nvmax = p, method = "backward")
    reg_summary_bwd = summary(out_bwd)
    
    reg_summary_bwd
    reg_summary_bwd$adjr2
    
    best_bwd = which.max(reg_summary_bwd$adjr2)
    
    best_bwd_pos = which(reg_summary$outmat[best_bwd,] == "*")
    
    gam_best_bwd = rep(0,p)
    gam_best_bwd[1] = 1
    
    gam_best_bwd[best_bwd_pos + 1] = 1
    gam_best_bwd
    
    out_best_bwd = lm(y ~.-1, data= X[,which(gam_best_bwd == 1)], subset = train)
    
    X_test = as.data.frame(X_test)
    
    y_hat_best_bwd = predict(out_best_bwd, newdata = X_test[,which(gam_best_bwd == 1)])
    
    mse_best_bwd = mean((y_test - y_hat_best_bwd)^2)
    MSE_bwd[i,] = mse_best_bwd
    
    ### LASSO ###
    
    grid = 10^seq(-3, 9, length = 100)
    
    X = as.matrix(X)
    cv.out = cv.glmnet(X, y, lambda = grid, alpha = 1)
    
    bestlam = cv.out$lambda.min
    
    lasso_best = glmnet(X, y, alpha = 1, lambda = bestlam, standardize = TRUE)
    lasso_pos = lasso_best$beta@i
    
    gam_lasso = rep(0, p)
    gam_lasso[1] = 1
    
    gam_lasso[lasso_pos + 1] = 1
    
    
    X = as.data.frame(X)
    
    out_lasso = lm(y ~.-1, data = X[,which(gam_lasso == 1), drop = FALSE], subset = train)
    
    X_test = as.data.frame(X_test)
    
    y_hat_lasso = predict(out_lasso, newdata = X_test[,which(gam_lasso == 1), drop = FALSE])
    
    mse_lasso = mean((y_test - y_hat_lasso)^2)
    MSE_lasso[i,] = mse_lasso
    
    ### FMA ###
    
    AICs = c()
    
    X = as.matrix(X)
    
    Y_pred = c()
    
    
    out_tmp = glm(y ~ .-1, data = as.data.frame(X[, 1, drop = FALSE]), subset = train)
    
    Y_pred = cbind(Y_pred, predict(out_tmp, newdata = X_test)) # use test inside predict
    
    t = out_tmp$aic
    
    names(t) = paste("0")
    
    AICs = c(AICs,t)
    
    ##
    
    for (q in 1:(p-1)) {
      all_q = as.matrix(combn(1:(p-1),q))
      
      for (k in 1:ncol(all_q)) {
        
        out_tmp = glm(y ~.-1, data = as.data.frame(X[,c(1,all_q[,k] + 1)]), subset = train) # use the train, di nuovo check problema nomi variabili
        
        Y_pred = cbind(Y_pred, predict(out_tmp, newdata = X_test)) # use test inside predict
        
        t = out_tmp$aic
        
        names(t) = paste(all_q[,k] + 1, collapse = ',')
        
        AICs = c(AICs,t)
      }
      
    }
    
    
    freq_weights = exp(-AICs/2)/sum(exp(-AICs/2))
    
    y_hat_fma = Y_pred%*%freq_weights
    
    mse_fma = mean((y_hat_fma - y_test)^2) # use test
    MSE_FMA[i,] = mse_fma
    
    X = NULL
    y = NULL
    
  }
  
  
  out_MSE = data.frame(bestsubset = MSE_best_subset,
            forward = MSE_forward,
            backward = MSE_bwd,
            lasso = MSE_lasso,
            FMA = MSE_FMA,
            BMA = MSE_BMA)
  
  return(out = list(out_MSE = out_MSE,
                    X_all = X_all))
  
}

### p = 5 ###

out_sim_n50_p5_pred = sim_over_samples_pred(n = 50, p = 5, N = 10)
out_sim_n100_p5_pred = sim_over_samples_pred(n = 100, p = 5, N = 10)
out_sim_n200_p5_pred = sim_over_samples_pred(n = 200, p = 5, N = 10)
out_sim_n500_p5_pred = sim_over_samples_pred(n = 500, p = 5, N = 10)
out_sim_n1000_p5_pred = sim_over_samples_pred(n = 1000, p = 5, N = 10)

### p = 10 ###

out_sim_n50_p10_pred = sim_over_samples_pred(n = 50, p = 10, N = 10)
out_sim_n100_p10_pred = sim_over_samples_pred(n = 100, p = 10, N = 10)
out_sim_n200_p10_pred = sim_over_samples_pred(n = 200, p = 10, N = 10)
out_sim_n500_p10_pred = sim_over_samples_pred(n = 500, p = 10, N = 10)
out_sim_n1000_p10_pred = sim_over_samples_pred(n = 1000, p = 10, N = 10)

### p = 15 ###

out_sim_n50_p15_pred = sim_over_samples_pred(n = 50, p = 15, N = 10)
out_sim_n100_p15_pred = sim_over_samples_pred(n = 100, p = 15, N = 10)
out_sim_n200_p15_pred = sim_over_samples_pred(n = 200, p = 15, N = 10)
out_sim_n500_p15_pred = sim_over_samples_pred(n = 500, p = 15, N = 10)
out_sim_n1000_p15_pred = sim_over_samples_pred(n = 1000, p = 15, N = 10)

# save.image("out_all_prediction.RData")

# Results

# p=5

round(colMeans(out_sim_n50_p5_pred$out_MSE), 2)
round(colMeans(out_sim_n100_p5_pred$out_MSE),2)
round(colMeans(out_sim_n200_p5_pred$out_MSE),2)
round(colMeans(out_sim_n500_p5_pred$out_MSE),2)
round(colMeans(out_sim_n1000_p5_pred$out_MSE),2)

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


#save.image("PredN.Rdata")


load("PredN.Rdata")

