###### estimate function ######
#gamma: parameter DPD function
#lambda : penalty function parameter 
#weigh.mode = Ad-DPD or AW-DPD
library(ncvreg)
library(robustHD)
library(glmnet)

dpd_adaptLASSO<- function(X, Y, beta, beta0, sigma, weight.mode = c("lasso","adapt","scad"), lambda, gamma, inter, regul_alpha){
  #if every coeff on init beta is zero, stop
  if(all(beta==0)){
    stop("null beta init- non penalized DPD")
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
    
    # #calculate weight
    w = rep(0,p)
    weight <- match.arg(weight.mode)
    
    if(weight =="lasso"){
      w = w+1/N
    }else if(weight == "adapt"){
      #calculate weights of adaptitive LASSO
      for(i in 1:p){if(beta[i] != 0){w[i] = 1/abs(beta[i])} }
      #for zero coefficients the greatest penalization
      w[w == 0] = 4*max(w)
    }else if(weight == "scad"){
      #else SCAD penalty
      a=3.7
      for(i in 1:p){
        if(abs(beta[i]) <= lambda){w[i] = 1
        }else if(lambda < abs(beta[i]) & abs(beta[i]) <= a*lambda){
          w[i] = (a*lambda - abs(beta[i]))/((a-1)*lambda)
        }
      }
    }else{stop("you need to specify a valid weight mode")}
    

    #temporary copies 
    beta0_tmp = beta0*inter
    beta_tmp = beta
    sigma_tmp = sigma
    
    #weights of MM-algorithm
    mu = exp(-(gamma/2)*((Y-tmp3)/sigma)^2)
    if(sum(mu) == 0){
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
    estimate <- glmnet(X_w, Y_w, intercept = FALSE, standardize = TRUE, alpha = regul_alpha,
                      family = "gaussian", lambda = lambda, penalty.factor = w)
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
                             + N*lambda*sum(abs(w*beta_tmp))- N*lambda*sum(abs(w*beta)) )
    if(all(beta==0) | stop2 < 1e-9 ){
      break
    }
  }
  return(list("beta0"=beta0,"beta" = beta, "sigma"= sigma))
}



