## generating outputs for multiple settings

p=500
n=100
epsilon = .1

intercept = 1
num_signal = 1

rm(list=ls())

setwd("G:/Mi unidad/UCM/ARTICULOS/Ad-DPD-LASSO/code/Codigos segunda version")
library(MASS)


library(parallel)
library(doParallel)

library(glmnet)
library(ncvreg)
library(robustHD)
library(rqPen)

library(tilting)
library(flare)

source("dpd_estimation.R")
source("cv_dpd_samelambda.R")



# function to get metrics
get.metrics = function(beta, beta.hat,  sigma, Xtest, Ytest){
  signal = which(beta!=0)
  #MS    TP   TN |  MSES   MSEN   EE | APrB
  c(sum(beta.hat != 0), # Model size
    sum(beta != 0 & beta.hat != 0)/sum(beta != 0), # True positive
    sum(beta == 0 & beta.hat == 0)/sum(beta == 0), # True negative
    mean((beta[signal]-beta.hat[signal])^2), # MSES = mean sq est error S
    mean((beta[-signal]-beta.hat[-signal])^2), # MSEN = mean sq est error N
    abs(.5-sigma),
    sqrt(mean(Ytest - Xtest %*% beta.hat)^2)) # BIAS
}


## Master simulation function
out.fun = function(n, p, num_signal, epsilon, nrep){
  
  ## set up simulation quantities
  tB = rep(0,p)
  # random_sample = sample(1:blocks,3)
  for(b in 1:num_signal){tB[(20*(b-1)+1):(20*(b-1)+5)] = c(3,1.5,0,0,2)}
  
  Sigma = diag(p)
  for(i in 2:p){
    for(j in 1:(i-1)){
      Sigma[i,j] = .5^abs(i-j)
      Sigma[j,i] = Sigma[i,j]
    }
  }
  
  ## function to run over repetitions
  loopfun = function(nrep){
    
    set.seed(nrep*1e3)
    X = mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
    Y = drop(X %*% tB + rnorm(n)*.5)
    
    ## add outliers on X
    nout = ceiling(100*epsilon)
    if (nout > 0){
      X[1:10,30:(29+nout)] = X[1:10,30:(29+nout)] + 20+rnorm(nout)}
      #Y[1:nout] = Y[1:nout] + 20 +rnorm(nout)}
    
    # test data
    Xtest = mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
    Ytest = drop(Xtest %*% tB + rnorm(n)*.5)
    
    err.mat1 = matrix(NA, ncol=8, nrow=32)
    
    gamma_list = c(0.1, 0.3,0.5, 0.7,1) 
    #gamma_list = c(0.5,0.7,1)
    for(i in 1:5){
      
      # if(i == 1 ){lmin = 3e-2; lmax= .09
      # }else if(i==2){lmin = 2.5e-2; lmax= .07
      # }else if(i==5){lmin = 1.9e-2; lmax= .06
      # }else{lmin = 2.1e-2; lmax=.09}
      
      if(i == 1 | i ==2 ){lmin = 2.7e-2; lmax= .08
      #}else if(i==2){lmin = 2.5e-2; lmax= .07
      #}else if(i==5){lmin = 1.9e-2; lmax= .06
      }else{lmin = 2.7e-2; lmax=.08}
      
      #LASSO
      t1 = system.time(z <- HBIC.dpd(X, as.matrix(Y), intercept=T,
                                       init.mode="RLARS", init.model=NULL, penalty = "adaptive",
                                       weight.mode = "lasso",lambda.mode="given",
                                       lmin=lmin, lmax=lmax, gam=gamma_list[i], nlambda = 50))

      if(!is.na(z$best$beta)){ err.mat1[i,] = c(get.metrics(tB, z$best$beta, z$best$sigma, Xtest, Ytest), t1[3]) }
      
      #ADAPTATIVE LASSO
      t2 = system.time(z <- HBIC.dpd(X, as.matrix(Y), intercept=T,
                                     init.mode="DPD-lasso", init.model=NULL, penalty = "adaptive",
                                     weight.mode = "adaptive",lambda.mode="lambda0",
                                     lmin=lmin, lmax=lmax, gam=gamma_list[i], nlambda = 50))

      if(!is.na(z$best$beta)){ err.mat1[5+i,] = c(get.metrics(tB, z$best$beta, z$best$sigma, Xtest, Ytest), t2[3]) }

      #ADAPTATIVE SCAD
      t3 = system.time(z <- HBIC.dpd(X, as.matrix(Y), intercept=T,
                                     init.mode="DPD-lasso", init.model=NULL, penalty = "adaptive",
                                     weight.mode = "scad",lambda.mode="given",
                                     lmin=lmin, lmax=lmax, gam=gamma_list[i], nlambda = 50))

      if(!is.na(z$best$beta)){ err.mat1[10+i,] = c(get.metrics(tB, z$best$beta, z$best$sigma, Xtest, Ytest), t3[3]) }

      #SCAD
      t4 = system.time(z <- HBIC.dpd(X, as.matrix(Y), intercept=T,
                                     init.mode="DPD-lasso", init.model=NULL, penalty = "SCAD",
                                     weight.mode = "lasso",lambda.mode="given",
                                     lmin=lmin, lmax=lmax, gam=gamma_list[i], nlambda = 50))

      if(!is.na(z$best$beta)){ err.mat1[15+i,] = c(get.metrics(tB, z$best$beta, z$best$sigma, Xtest, Ytest), t4[3]) }

      #MCP
      t5 = system.time(z <- HBIC.dpd(X, as.matrix(Y), intercept=T,
                                     init.mode="DPD-lasso", init.model=NULL, penalty = "MCP",
                                     weight.mode = "lasso",lambda.mode="given",
                                     lmin=lmin, lmax=lmax, gam=gamma_list[i], nlambda = 50))

      if(!is.na(z$best$beta)){ err.mat1[20+i,] = c(get.metrics(tB, z$best$beta, z$best$sigma, Xtest, Ytest), t5[3]) }

    }
    
    ## LS-LASSO
    t1 = system.time(z <- cv.glmnet(X, as.matrix(Y), nfolds=5))
    err.mat1[26,] = c(get.metrics(tB, coef(z)[-1], 1.4826*mad(Y - X %*% coef(z)[-1]), Xtest, Ytest), t1[3])
    
    ## LS-SCAD
    t2 = system.time(z <- cv.ncvreg(X, Y, family="gaussian", penalty="SCAD"))
    beta.SCAD = z$fit$beta[,which.min(z$cve)][-1]
    err.mat1[27,] = c(get.metrics(tB, beta.SCAD, 1.4826*mad(Y - X %*% beta.SCAD), Xtest, Ytest), t2[3])

    ## LS-MCP
    t3 = system.time(z <- cv.ncvreg(X, Y, family="gaussian", penalty="MCP"))
    beta.MCP = z$fit$beta[,which.min(z$cve)][-1]
    err.mat1[28,] = c(get.metrics(tB, beta.MCP, 1.4826*mad(Y - X %*% beta.MCP), Xtest, Ytest), t3[3])

    # ## LAD-lasso
    # t4 = system.time(LAD <- slim(X, Y, nlambda=50, q=1, lambda.min.value=.01, verbose = FALSE))
    # BICvals = as.numeric(with(LAD, lapply(1:nlambda,
    #                                       function(i) sum(abs(Y - X %*% beta[,i])) +
    #                                         log(n)*sum(beta[,i]!=0))))
    # beta.LAD = LAD$beta[,which.min(BICvals)]
    # err.mat1[29,] = c(get.metrics(tB, beta.LAD, 1.4826*mad(Y - X %*% beta.LAD), Xtest, Ytest), t4[3])
    
    #rm(LAD)
    # RLARS
    t5 = system.time(RLARS <- rlars(Y~X))
    err.mat1[30,] = c(get.metrics(tB, coef(RLARS)[-1], getScale(RLARS), Xtest, Ytest), t5[3])
    rm(RLARS)
    
    # # sLTS
    # frac<- seq(0,lambda0(X0,Y) ,by=(1/50)*lambda0(X0,Y))
    # t26 = system.time(sLTS <- sparseLTS(Y~X0, inte=FALSE, lambda=frac[-1], mode="lambda"))
    # err.mat1[31,] = c(get.metrics(tB, coef(sLTS)[-c(1,2)], getScale(sLTS), Xtest, Ytest), t26[3])
    # 
    # rm(sLTS)
    # 
    # # RANSAC
    # t47 = system.time(ransac <- gam.ini(X0, as.matrix(Y), ini.subsamp=0.2, cl.num=1, ini.cand=1000, reg.alpha=1))
    # err.ma1t[32,] = c(get.metrics(tB, as.matrix(ransac$beta[-1]), ransac$sigma, Xtest, Ytest), t47[3])
    # rm(ransac)
    gamma_list = c(0.1, 0.3,0.5, 0.7,1) 
    colnames(err.mat1) <- c("Size","TP","TN","MSES", "MSEN", "EE(Sigma)","APrB","Runtime")
    row.names(err.mat1) <- c(paste0("DPD-LASSO gamma=", gamma_list), 
                              paste0("Ad-DPD-LASSO gamma=", gamma_list),
                              paste0("AW-DPD-LASSO gamma=", gamma_list),
                              paste0("DPD-SCAD gamma=", gamma_list),
                              paste0("DPD-MCP gamma=", gamma_list),
                              "LS-LASSO","LS-SCAD","LS-MCP","LAD-Lasso",
                              "RLARS","sLTS","RANSAC")
    return(err.mat1) } 
  
  # # parallel programming
  numCores <- 3
  cl <- makeCluster(numCores)
  clusterExport(cl,c("get.metrics"))
  clusterExport(cl, c("n","p","epsilon"), envir= environment())
  clusterEvalQ(cl,c("tB","Sigma"))
  clusterCall(cl = cl,function() library(MASS))
  clusterCall(cl = cl,function() library(robustHD))
  clusterCall(cl = cl,function() library(ncvreg))
  clusterCall(cl = cl,function() library(rqPen))
  
  clusterCall(cl = cl,function() library(glmnet))
  # clusterCall(cl = cl,function() library(tilting))
  # clusterCall(cl = cl,function() library(flare))
  
  clusterCall(cl, function() {source("dpd_estimation.R") })
  clusterCall(cl, function() {source("cv_dpd_samelambda.R") })
  
  out.list<- parLapply(cl,1:nrep,loopfun)
  stopCluster(cl)
  #out.list = lapply(1:nrep,loopfun)
  return(out.list)}

mean_error <- function(out.list){
  
  err = data.frame(matrix(0, ncol=8, nrow=32))
  for (i in 1:length(out.list)){err = err + out.list[[i]] }
  err = err/length(out.list)
  names(err) = c("Size","TP","TN","MSES", "MSEN", "EE(Sigma)","APrB","Runtime")
  gamma_list = c(0.1, 0.3,0.5, 0.7,1) 
  row.names(err) = c(paste0("DPD-LASSO gamma=", gamma_list), 
                     paste0("Ad-DPD-LASSO gamma=", gamma_list),
                     paste0("AW-DPD-LASSO gamma=", gamma_list),
                     paste0("DPD-SCAD gamma=", gamma_list),
                     paste0("DPD-MCP gamma=", gamma_list),
                     "LS-LASSO","LS-SCAD","LS-MCP","LAD-Lasso",
                     "RLARS","sLTS","RANSAC" )
  return(err)
}

out500e0n1_list = out.fun(n=100, p=500, num_signal = 1, epsilon = 0, nrep=100)
out500e0n1 = mean_error(out500e0n1_list)

out500e0n3_list = out.fun(n=100, p=500, num_signal = 3, epsilon = 0, nrep=100)
out500e0n3 = mean_error(out500e0n3_list)

#out
out500e1n1_list = out.fun(n=100, p=500, num_signal = 1, epsilon = .1, nrep=20)
out500e1n1 = mean_error(out500e1n1_list)
# 
out500e1n3_list = out.fun(n=100, p=500, num_signal = 3, epsilon = .1, nrep=20)
out500e1n3 = mean_error(out500e1n3_list)

# #1000
out1000e0n1_list = out.fun(n=100, p=1000, num_signal = 1, epsilon = 0, nrep=100)
out1000e0n1 = mean_error(out1000e0n1_list)
# 
out1000e0n3_list = out.fun(n=100, p=1000,  num_signal = 3, epsilon = 0, nrep=100)
out1000e0n3 = mean_error(out1000e0n3_list)

#out
out1000e1n1_list = out.fun(n=100, p=1000, num_signal = 1, epsilon = .1, nrep=100)
out1000e1n1 = mean_error(out1000e1n1_list)

out1000e1n3_list = out.fun(n=100, p=1000,  num_signal = 3, epsilon = .1, nrep=100)
out1000e1n3 = mean_error(out1000e1n3_list)

save(out500e0n1, file = "outDPDp500e0n1.Rda")
save(out500e0n3, file = "outDPDp500e0n3.Rda")
save(out500e1n1, file = "outDPDp500ex1n1.Rda")
save(out500e1n3, file = "outscadDPDp500ex1n3.Rda")

save(out1000e0n1, file = "outDPDp1000e0n1.Rda")
save(out1000e0n3, file = "outDPDp1000e0n3.Rda")
save(out1000e1n1, file = "outDPDp1000ex1n1.Rda")
save(out1000e1n3, file = "outDPDp1000ex1n3.Rda")


