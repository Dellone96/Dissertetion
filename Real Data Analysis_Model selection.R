library(invgamma)
library(mvtnorm)

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

### BMA ###

#### Marginal Distribution ####

beta0  = rep(0, p)
a = 0.01
b = 0.01
V = diag(1,p)

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


final = c()


for (q in 1:10) { # q è il numero di variabili inserite nel modello, va da 1 fino a un massimo di 10
  all_q = combn(1:p,q) # qui trovo le combinazioni delle p variabili (1:p) prese q alla volta
  for (n in 1:ncol(all_q)) {
    t = m_y(y = y, X = X[,c(1,all_q[,n] + 1)], a = a, b = b, beta0 = beta0, V = V)  
    names(t) = paste(all_q[,n] + 1, collapse = ',')
    
    final = c(final,t)
  }
  print(q)
}

final

post = exp(final)/sum(exp(final))

c = mean(final)
c

post_c = exp(final - c)/(sum(exp(final - c)))


post
post_c


#### Best bayesian model ####

## This is the Maximum A Posteriori Model

best_bayes = which.max(post_c)

post_c[best_bayes]

## Compute posterior probabilities of inclusion, for each variable in the dataset

## And construct the model which includes those predictors whose posterior probability is > 0.5
## This is called Median Probability Model (MPM)

tab_models = matrix(0, nrow = length(post_c), ncol = p+1)

tab_models[,1] = 1

for(i in 1:length(post_c)){
  
  tab_models[i, as.numeric(unlist(strsplit(x = names(post_c[i]), split = ",")))] = 1

}

probs_inclusion = c(round(t(tab_models)%*%post_c, 3))
probs_inclusion

names(probs_inclusion) = c(colnames(X)[-1])

par(mar = c(7,3,1,1))

barplot(c(probs_inclusion), las = 2)


gamma_bayes = (probs_inclusion > 0.5)*1

gamma_bayes


#########
#########


sum(post)
barplot(post)

summary(lm(y ~ X - 1))

gam_bayes = rep(0,p)

gam_bayes[as.integer(unlist(strsplit (names(best_bayes), ",")))]= 1

gam_bayes[1] = 1

#### Best Subset ####

library(ISLR)
library(leaps)

X = as.data.frame(X)


out_subset  = regsubsets(y ~ ., X[,-1], nvmax = p)
reg_summary = summary(out_subset)
reg_summary

reg_summary$adjr2
reg_summary$rss


best_subset = which.max(reg_summary$adjr2)

best_subset_pos = which(reg_summary$outmat[best_subset,] == "*")

gam_best_subset = rep(0,p)
gam_best_subset[1] = 1

gam_best_subset[best_subset_pos + 1] = 1
gam_best_subset

#### Forward Stepwise ####

out_fwd = regsubsets(y ~ ., X, nvmax = p, method = "forward")
reg_summary_fwd = summary(out_fwd)

reg_summary_fwd
reg_summary_fwd$adjr2

best_forward = which.max(reg_summary_fwd$adjr2)

best_forward_pos = which(reg_summary$outmat[best_forward,] == "*")

gam_best_forward = rep(0,p)
gam_best_forward[1] = 1

gam_best_forward[best_forward_pos + 1] = 1
gam_best_forward

#### Backward ####

out_bwd = regsubsets(y ~ ., X, nvmax = p, method = "backward")
reg_summary_bwd = summary(out_bwd)

reg_summary_bwd
reg_summary_bwd$adjr2

best_bwd = which.max(reg_summary_bwd$adjr2)

best_bwd_pos = which(reg_summary$outmat[best_bwd,] == "*")

gam_best_bwd = rep(0,p)
gam_best_bwd[1] = 1

gam_best_bwd[best_bwd_pos + 1] = 1
gam_best_bwd

#### Lasso ####

library(glmnet)

set.seed(1)

#grid value of lambda

grid = 10^seq(-3, 9, length = 100)

X = as.matrix(X)
cv.out = cv.glmnet(X, y, lambda = grid, alpha = 1)


bestlam = cv.out$lambda.min
bestlam

lasso_best = glmnet(X, y, alpha = 1, lambda = bestlam, standardize = TRUE)
lasso_pos = lasso_best$beta@i

gam_lasso = rep(0, p)
gam_lasso[1] = 1

gam_lasso[lasso_pos + 1] = 1

### FMA ###

AICs = c()

for (q in 1:10) { # q è il numero di variabili inserite nel modello, va da 1 fino a un massimo di 10
  all_q = combn(1:p,q) # qui trovo le combinazioni delle p variabili (1:p) prese q alla volta
  for (k in 1:ncol(all_q)) {
    
    out_tmp = glm(y ~ X[,c(1,all_q[,k] + 1)] - 1)
    
    t = out_tmp$aic
    
    names(t) = paste(all_q[,k] + 1, collapse = ',')
    
    AICs = c(AICs,t)
  }
  print(q)
}


freq_weights = exp(-AICs/2)/sum(exp(-AICs/2))


best_FMA = which.min(AICs)

gam_FMA = rep(0,p)

gam_FMA[as.integer(unlist(strsplit (names(best_FMA), ",")))]= 1

gam_FMA[1] = 1


## Posterior probabilities of inclusion (BMA)

tab_models = matrix(0, nrow = length(post_c), ncol = p+1)

tab_models[,1] = 1

for(i in 1:length(post_c)){
  
  tab_models[i, as.numeric(unlist(strsplit(x = names(post_c[i]), split = ",")))] = 1
  
}

probs_inclusion = c(round(t(tab_models)%*%post_c, 3))
probs_inclusion

colnames(X)[1] = "Intercept"

names(probs_inclusion) = c(colnames(X))

par(mar = c(7,3,1,1))

barplot(c(probs_inclusion), las = 2)

gamma_bayes = (probs_inclusion > 0.5)*1
gamma_bayes


## Posterior probabilities of inclusion (FMA)

tab_models = matrix(0, nrow = length(freq_weights), ncol = p+1)

tab_models[,1] = 1

for(i in 1:length(freq_weights)){
  
  tab_models[i, as.numeric(unlist(strsplit(x = names(freq_weights[i]), split = ",")))] = 1
  
}

probs_inclusion_fma = c(round(t(tab_models)%*%freq_weights, 3))
probs_inclusion_fma

names(probs_inclusion_fma) = c(colnames(X))

par(mar = c(7,3,1,1))

barplot(c(probs_inclusion_fma), las = 2)

gamma_fma = (probs_inclusion_fma > 0.5)*1
gamma_fma

#save.image("RD models.Rdata")

load("RD models.Rdata")