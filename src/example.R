rm(list=ls())

###################################
# example for logistic regression #
###################################

library(MASS)
library(glmnet)

source("dpd_estimation.R")
source("logit.R")
source("hgic_dpd_logit.R")

k=500 #number of variables
n=100 #sample size
alpha = 0.3 #DPD loss tunning parameter

num_signal = 1 #number of nonzero coefficients
epsilon = 0.05 #outliers percent


## True regression parameters
tB = rep(0,k)
for(b in 1:num_signal){tB[(20*(b-1)+1):(20*(b-1)+5)] = c(5,5,0,0,5)}

## Variance-covariance matrix of explanatory variables
Sigma = diag(k)
  for(i in 2:k){
    for(j in 1:(i-1)){
      Sigma[i,j] = 0.1^abs(i-j)
      Sigma[j,i] = Sigma[i,j]
    }
  }
  
#simulate data
X = mvrnorm(n, mu=rep(0,k), Sigma=Sigma)

pi = sapply(drop(X %*% tB), invlogit)
pi[pi<1e-4]=0; pi[pi>(1-1e-4)]=1 ###avoiding large estimates 
    
if(epsilon >0){
  raffle = rbinom(size=1,n=n,p=(1-epsilon))
  Y1 =sapply(pi,function(x) rbinom(size=1,n=1,p=x)) # pure data
  Y2 = sapply(pi,function(x) rbinom(size=1,n=1,p=(1-x))) # contamination
  Y  = raffle*Y1+(1-raffle)*Y2
}else{Y = sapply(pi,function(x) rbinom(size=1,n=1,p=x))}

X=cbind(1,X)

# initial estimator
LS_LASSO<-cv.glmnet(X[,2:(k+1)], Y, family = "binomial", intercept=T)
beta.init = coef(LS_LASSO)[1:(k+1)]
beta.init[beta.init!=0]

#different adaptive estimators based on DPD loss 
lmin = 0.005;lmax=0.5
#lasso
z <- HGIC.dpdLASSO_LR(X,Y, beta.init, weight.mode = "lasso",LR.loss = "LS", lambda.mode="given",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]

#adaptive
z <- HGIC.dpdLASSO_LR(X,Y, beta.init, weight.mode = "adaptive", LR.loss = "LS", lambda.mode="lambda0",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]
      
#scad
z <- HGIC.dpdLASSO_LR(X,Y, beta.init, weight.mode = "scad", LR.loss = "LS", lambda.mode="lambda0",lmax=lmax,lmin=lmin, nlambda=50, alpha =alpha)
z$best$beta[z$best$beta !=0]