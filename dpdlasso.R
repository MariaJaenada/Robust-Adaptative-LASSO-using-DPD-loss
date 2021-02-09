###### estimate function ######
#gamma: parameter Renyi function
#lambda : penalty function parameter 
library(ncvreg)
library(robustHD)
library(glmnet)

dpd_LASSO<- function(X, Y, beta, beta0, sigma, lambda, gamma, inter, regul_alpha){
  #if every coeff on init beta is zero, stop
  if(all(beta==0)){
    warning("null beta init- non penalized DPD")
  }
  
  N = dim(X)[1]
  p = dim(X)[2]
  
  #create intercept term
  tmp = rep(1, N)
  tmp1 = inter*beta0*tmp
  tmp2 = X%*%beta #X matrix doesnt contain ones column 
  tmp3 = drop(tmp1 + tmp2) #y estimate
  
  # iterative minimization 
  for (m in 1:500){
    #temporary copies 
    beta0_tmp = beta0*inter
    beta_tmp = beta
    sigma_tmp = sigma
    
    #weights of MM-algorithm
    mu = exp(-(gamma/2)*((Y-tmp3)/sigma)^2)
    if(sum(mu)==0){
      beta0 = "NaN"
      beta = "NaN"
      sigma = "NaN"
      break}
    mu = drop(mu/sum(mu))
    
    #update beta0: doesn't depend on penalization
    beta0 = drop(t(mu)%*%(Y - tmp2)*inter)
    tmp1 = beta0*tmp
    
    #weigthed matrices : reduce to least squares problem
    Y_w = diag(sqrt(mu))%*%((Y-tmp1)/sigma) #Y contains incercept (mean)
    X_w =  diag(sqrt(mu))%*%(X/sigma) #without intercept
    
    #minimize beta with adaptative lasso
    #Solve using LARS
    estimate <- glmnet(X_w, Y_w, intercept = FALSE, standardize = TRUE, alpha = regul_alpha,
                       family = "gaussian", lambda = lambda)
    
    beta = coef(estimate)[1:p+1]
  
    #update tmp3 with new beta
    tmp2 = X%*%beta 
    tmp3 = drop(tmp1 + tmp2) 
    
    #update sigma now
    wi = exp((-gamma/2)*((Y-tmp3)/sigma)^2)
    sigma = sqrt(abs(drop(mean(wi*(Y-tmp3)^2))*(mean(wi)-gamma/(gamma+1)^(3/2))^(-1)))
    
    if(sigma < 1e-7){
      beta0 = "NaN"
      beta = "NaN"
      sigma = "NaN"
      break}

    ##stopping criteria LASSO penalization
    stop2 = abs(sigma^(-gamma)*(1/(1+gamma)^(3/2)-(1/(N*gamma)*sum(exp(-(gamma/2)*((Y-tmp3)/sigma)^2)) ))
                - (sigma_tmp^(-gamma)*(1/(1+gamma)^(3/2)-(1/(N*gamma)*sum(exp(-(gamma/2)*((Y-beta0_tmp*tmp-X%*%beta_tmp)/sigma_tmp)^2)) )))
                + N*lambda*sum(abs(beta_tmp))- N*lambda*sum(abs(beta)) )

    if(all(beta==0) | stop2 < 1e-9 ){
      break
    }
  }
  return(list("beta0"=beta0,"beta" = beta, "sigma"= sigma))
}


HBIC.dpdLASSO <- function(X,Y,lambda.mode = "given",
                          lmax=1,lmin=0.05, nlambda=50,
                          ncores=1, gam=0.1, gam0=0.5, intercept="TRUE", alpha=1,
                          ini.subsamp=0.2, ini.cand=1000,
                          alpha.LTS=0.75, nlambda.LTS=40){
  
  
  #########################################
  # detect the data type and tuning parameters value
  #########################################
  if(is.matrix(X)=="FALSE" || is.matrix(Y)=="FALSE" ){
    stop(message="The data type of X and Y is only Matrix.")
  }
  if(gam <0 || gam0 <0){
    stop(message="Robust tuning parameters gam and gam0 are only positive.")
  }

  ini.cand <- as.integer(ini.cand)
  if(ini.cand <2  ){
    stop(message="The number of candidates of initial points is more than 2.")
  }
  nlambda <- as.integer(nlambda)
  if(nlambda < 1 ){
    stop(message="The number of grids for sparse tuning parameter is more than 1.")
  }
 
  if(lmin < 0 || lmax <0){
    stop(message="lmin and lmax are positive.")
  }
  ncores <- as.integer(ncores)
  if(ncores < 1){
    stop(message="The number of cpu cores for parallel computing is more than 1.")
  }
  if( ini.subsamp <0 || ini.subsamp >1){
    stop(message="The proportion of subsample is 0 < ini.subsamp < 1")
  }
  
  #########################################
  # intercept
  #########################################
  
  if(intercept=="TRUE"){
    inter <-1
  }else {
    inter <-0
  }
  
  #########################################
  # choose a initial point
  ########################################
  
  init.time = system.time(RLARS <- rlars(Y~X, seed = 1024))
  beta0.init <- inter*(coef(RLARS)[1])
  beta.init <- as.matrix(coef(RLARS)[-1])
  sigma.init <- getScale(RLARS)
  sigma.cv <- sigma.init
 
  
  ########################################
  # regularization grids
  ########################################
  
  if(lambda.mode=="lambda0"){
    tmp <- lambda.ini(X, Y, beta0.init, beta.init, sigma.init ,gam)
    lambda <- exp(seq(  log(0.05*tmp)  ,log(tmp)  ,length=nlambda))
  }else{
    lambda <- exp(seq( log(lmin)  ,log(lmax)  ,length=nlambda))
  }
  
  #######################################
  # caluculate
  #######################################
  n = nrow(X); p = ncol(X)
  Cn = log(log(n))
  mult = Cn*log(p)/n
  
  
  # calculate solutions
  res.reggam <- list()
  HBIC.vec = rep("NaN", nlambda)
  for(k in 1:nlambda){
    
    res <- dpd_LASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init,
                          lambda[k], gam, inter, 1)
    if(res$beta0 != "NaN"){
      res.reggam[[k]] <- list(beta0=res$beta0, beta=res$beta, sigma=res$sigma)
      M = sum(res$beta != 0)
      SSE = res$sigma^2 # robust measure of variance
      HBIC.vec[k] = log(SSE) + M*mult
    }
  }
  
  
  # for large gamma, some values in HBIC.vec are NaN so control for that
  best=NA # if all are NaN then no best model
  which.good = which(paste(HBIC.vec)!="NaN")
  if(length(which.good)>0){
    best_HBIC = min(HBIC.vec[which.good])
    best = res.reggam[[which(HBIC.vec == best_HBIC)[1]]] # get good entries and best model among them
  }
  print(lambda[[which(HBIC.vec == best_HBIC)[1]]])
  
  
  return(list(init.time=init.time[3], lambda=lambda, error=HBIC.vec,
              fits=res.reggam, best=best))
  
}
