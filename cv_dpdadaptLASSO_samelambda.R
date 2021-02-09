HBIC.dpdadaptLASSO <- function(X,Y,
                               init.mode=c("RLARS","DPD-lasso"), init.model=NULL,
                               weight.mode = "lasso",
                               lambda.mode=c("lambda0","given"),
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
  #########################################
  init <- match.arg(init.mode)
  
  if(init =="RLARS"){
    RLARS = init.model
    init.time = NULL
    if(is.null(RLARS)){
      init.time = system.time(RLARS <- rlars(Y~X, seed=1024))}
    
    beta0.init <- inter*(coef(RLARS)[1])
    beta.init <- as.matrix(coef(RLARS)[-1])
    sigma.init <- getScale(RLARS)
    sigma.cv <- sigma.init
  }else if(init == "DPD-lasso"){
    RLARS = init.model
    if(is.null(RLARS)){init.time = system.time(RLARS <- rlars(Y~X,seed = 1024))}
    beta0.init.lasso <- inter*(coef(RLARS)[1])
    beta.init.lasso <- as.matrix(coef(RLARS)[-1])
    sigma.init.lasso <- getScale(RLARS)
     
   }else{
    X1 = cbind(1,X)
    b = ginv(crossprod(X1)) %*% t(X1) %*% Y
    beta0.init = inter*b[1]
    beta.init = as.matrix(b[-1])
    sigma.init = 1.4826*mad(Y - X1 %*% b)
  }
  
  ########################################
  # regularization grids
  ########################################
  lmode <- match.arg(lambda.mode)
  
  if(lmode=="lambda0"){
    lambdamax0<-lambda0(X,Y) #Initial candidate
    lambdamax1<-try(optim.lam(X,Y,lambdamax0, gam, weight.mode))#Improvement
    if(class(lambdamax1)=='try-error'){
      lambdamax1<-lambdamax0
      warning('Using approximate lambdamax')
    }
    lambda<-seq(1e-8,lambdamax1,length.out=nlambda)
    if(dim(X)[2]>=dim(X)[1]){lambda<-lambda[2:nlambda]}
    
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
  for(k in 1:(nlambda-1)){
    #initialize with the same lambda
    if(init == "DPD-lasso"){
      DPD <- dpd_LASSO(X, as.matrix(Y), beta.init.lasso, beta0.init.lasso,
                     sigma.init.lasso, lambda[k], gam, inter, 1)
      beta.init <- DPD$beta
      beta0.init <- DPD$beta0
      sigma.init <- DPD$sigma
      if(beta0.init != "NaN"){
        if(sum(beta.init) !=0 ){
        res <- dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode = weight.mode,
                          lambda[k], gam, inter, 1)
        }else{res<- list(beta0 = "NaN")}
      }else{res <- list(beta0 = "NaN")}
    }else{
      res <- dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode = weight.mode,
                            lambda[k], gam, inter, 1)
    }
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

res.cal <- function(x,a){
  h <- floor(length(x)*a)
  t <- sqrt(sum(sort(x^2)[1:h])/h)
  return(t)
}
gam.ini <- function(X, Y, ini.subsamp, cl.num, ini.cand, reg.alpha){
  
  cl.ini <- makeCluster(cl.num)
  registerDoParallel(cl.ini)
  on.exit(stopCluster(cl.ini))
  ini.sam <- floor(nrow(X)*ini.subsamp)
  
  
  p <- ncol(X)
  n <- nrow(X)
  tmp <- foreach(k =1:ini.cand, .combine=rbind, .packages='glmnet', .export="res.cal")%dopar%{
    
    t <- sample(n,ini.sam)
    x <- X[t,]
    y <- Y[t,]
    
    x_t <- X[setdiff(1:n,t),]
    y_t <- Y[setdiff(1:n,t),]
    
    res <- glmnet(x,y,family="gaussian",alpha=reg.alpha)
    sigma <- apply( abs(matrix(y,nrow =nrow(x),ncol=length(res$lambda))- matrix(res$a0,nrow =nrow(x),ncol=length(res$lambda),byrow=TRUE)-x%*%res$beta ),2, median )
    sigma <- 1.4826*sigma
    tmp1 <- apply( matrix(y_t,nrow =nrow(x_t),ncol=length(res$lambda))- matrix(res$a0,nrow =nrow(x_t),ncol=length(res$lambda),byrow=TRUE)-x_t%*%res$beta,2, res.cal, 0.5 )
    tmp2 <-  which.min(tmp1)
    return(  c(tmp1[tmp2],res$a0[tmp2],res$beta[,tmp2], sigma[tmp2]) )
  }
  
  ini.ind <-  which.min(tmp[,1])
  return( list(beta0=unname(tmp[ini.ind, 2]), beta=unname(tmp[ini.ind, 3:(p+2)]), sigma=unname(tmp[ini.ind, (p+3)]))     )
}


optim.lam<-function(X,Y,lambdamax, gam, weight.mode){
  #Finds (approximately) smallest lambda such that the slope estimate is zero
  #INPUT
  #X, y: data set
  #lambdamax: initial guess for lambda
  #gam : parameter of DPD loss
  #weight.mode : scad or adapt
  #OUTPUT
  #lambda: smallest lambda such that the slope estimate is zero
  
  n<-nrow(X)
  p<-ncol(X)
  inter <- 1
  lambda2<-lambdamax
  RLARS <- rlars(Y~X, seed = 1024)
  
  beta0.init <- inter*(coef(RLARS)[1])
  beta.init <- as.matrix(coef(RLARS)[-1])
  sigma.init <- getScale(RLARS)
  
  beta.n<- dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode = weight.mode,
                          lambda2, gam, inter, 1)$beta
  zeros.n<-sum(beta.n==0)
  while (zeros.n<p & lambda2<max(n,1e4)){
    lambda2<-2*lambda2
    beta.n<-dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode = weight.mode,
                           lambda2, gam, inter, 1)$beta
    zeros.n<-sum(beta.n==0)
  }
  if(lambda2>=max(n,1e4))
  {
    return(lambda2)
  }
  lambda1<-lambdamax/2
  beta.o<-dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode =  weight.mode,
                         lambda1, gam, inter, 1)$beta
  zeros.o<-sum(beta.o==0)
  while(zeros.o==p){
    lambda1<-lambda1/2
    beta.o<-dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode =  weight.mode,
                           lambda1, gam, inter, 1)$beta
    zeros.o<-sum(beta.o==0)
  }
  for (j in 1:5){
    lambda<-0.5*(lambda1+lambda2)
    beta.n<-dpd_adaptLASSO(X, as.matrix(Y), beta.init, beta0.init, sigma.init, weight.mode =  weight.mode,
                           lambda, gam, inter, 1)$beta
    zeros.n<-sum(beta.n==0)
    if (zeros.n<p){
      lambda1<-lambda
    }else{
      lambda2<-lambda
    }
  }
  return(lambda2)
}


