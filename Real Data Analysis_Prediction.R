library(BMS)

data = as.data.frame(datafls)
dim(datafls)
n = nrow(data)
n

X = data[,c(6,11:16,25,28:35,38:40,42)]
X = X[,-c(6,10,12,13)]
y = data[,1]
dim(X)
p = ncol(X)


X = as.matrix(cbind(1,X))

head(X)

set.seed(123)

train = sample(1:n, n/2) # these observations will be used to fit the model
test = setdiff(1:n, train) # difference between sets 1:n and train

X_train = X[train,]
X_test  = X[test,]

y_train = y[train]
y_test = y[test]

#Best Subset
library(ISLR)
library(leaps)

X = as.data.frame(X)

out_subset_RD  = regsubsets(y ~ .-1, X, nvmax = p)

reg_summary_RD = summary(out_subset_RD)
reg_summary_RD

reg_summary_RD$adjr2

best_subset_RD = which.max(reg_summary_RD$adjr2)

gam_best_subset = rep(0,p)
gam_best_subset[1] = 1

best_subset_pos = which(reg_summary_RD$outmat[best_subset_RD,] == "*")

gam_best_subset[best_subset_pos + 1] = 1
gam_best_subset

## Evaluate the accuracy by fitting the model on the train and then using the test

out_best_subset = lm(y ~.-1, data= X[,which(gam_best_subset == 1)], subset = train)

X_test = as.data.frame(X_test)

y_hat_best_subset = predict(out_best_subset, newdata = X_test[,which(gam_best_subset == 1)])

mse_best_subset_RD = mean((y_test - y_hat_best_subset)^2)
mse_best_subset_RD*100

#### Forward Stepwise ####

out_fwd_RD = regsubsets(y ~ .-1, X, nvmax = p, method = "forward")
reg_summary_fwd_RD = summary(out_fwd_RD)

reg_summary_fwd_RD
reg_summary_fwd_RD$adjr2

best_forward_RD = which.max(reg_summary_fwd_RD$adjr2)

best_forward_pos = which(reg_summary_fwd_RD$outmat[best_forward_RD,] == "*")

gam_best_forward = rep(0,p)
gam_best_forward[1] = 1

gam_best_forward[best_forward_pos + 1] = 1
gam_best_forward

out_best_forward_RD = lm(y ~.-1, data= X[,which(gam_best_forward == 1)], subset = train)

X_test = as.data.frame(X_test)

y_hat_best_forward_RD = predict(out_best_forward_RD, newdata = X_test[,which(gam_best_forward == 1)])

mse_best_forward_RD = mean((y_test - y_hat_best_subset)^2)
mse_best_forward_RD*100

#### Backward ####

out_bwd_RD = regsubsets(y ~ .-1, X, nvmax = p, method = "backward")
reg_summary_bwd_RD = summary(out_bwd_RD)

reg_summary_bwd_RD
reg_summary_bwd_RD$adjr2

best_bwd_RD = which.max(reg_summary_bwd_RD$adjr2)

best_bwd_pos = which(reg_summary_RD$outmat[best_bwd_RD,] == "*")

gam_best_bwd = rep(0,p)
gam_best_bwd[1] = 1

gam_best_bwd[best_bwd_pos + 1] = 1
gam_best_bwd

out_best_bwd_RD = lm(y ~.-1, data= X[,which(gam_best_bwd == 1)], subset = train)

X_test = as.data.frame(X_test)

y_hat_best_bwd_RD = predict(out_best_bwd_RD, newdata = X_test[,which(gam_best_bwd == 1)])

mse_best_bwd_RD = mean((y_test - y_hat_best_bwd_RD)^2)
mse_best_bwd_RD*100

#### Lasso ####

library(glmnet)

#grid value of lambda

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

out_lasso_RD = lm(y ~.-1, data= X[,which(gam_lasso == 1)], subset = train)

X_test = as.data.frame(X_test)

y_hat_lasso_RD = predict(out_lasso_RD, newdata = X_test[,which(gam_lasso == 1)])

mse_lasso_RD = mean((y_test - y_hat_lasso_RD)^2)
mse_lasso_RD*100


# FMA Prediction

Y_pred = c()
AICs_pred = c()

for (q in 1:(p-1)) {
  all_q = as.matrix(combn(1:(p-1),q))
  
  for (k in 1:ncol(all_q)) {
    
    
    out_tmp = glm(y ~.-1, data = as.data.frame(X[,c(1,all_q[,k] + 1)]), subset = train) # use the train, di nuovo check problema nomi variabili
    
    Y_pred = cbind(Y_pred, predict(out_tmp, newdata = X_test)) # use test inside predict
    
    t = out_tmp$aic
    
    names(t) = paste(all_q[,k] + 1, collapse = ',')
    
    AICs_pred = c(AICs_pred,t)
  }
  
  print(q)
  
}


freq_weights = exp(-AICs_pred/2)/sum(exp(-AICs_pred/2))

y_hat_fma = Y_pred%*%freq_weights

mse_fma = mean((y_hat_fma - y_test)^2) # use test
mse_fma*100

#Bayesian model selection

source("functions_for_bma_prediction.R")

p = ncol(X)

beta0 = beta.0 = rep(0, p)
a = 0.01
b = 0.01
V = diag(1,p)

S = 500

out_k = posterior_prediction_linear_conjugate(y_train, y_test, X_train, X_test, beta.0, V, a, b, S)

str(out_k)

rowMeans(out_k$Y_pred)


## Repeat for all possible models corresponding to subsets of columns of X
## that is X[,c(1,all_q[,k] + 1)] inside the loop (as in the FMA with glm)

# combine all the model 


Y_pred_bma = c()

final = c()


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
  print(q)
}

post = exp(final)/sum(exp(final))

c = mean(final)

post_c = exp(final - c)/(sum(exp(final - c)))

post_c # weights

y_hat_bma = Y_pred_bma%*%post_c

mse_bma = mean((y_hat_bma - y_test)^2)
mse_bma*100

#save.image("RD pred.Rdata")

load("RD pred.Rdata")
