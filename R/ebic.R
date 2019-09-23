ebic_fl = function(x, y, beta, sigma2){
  n = nrow(x)
  p = ncol(x)
  beta <- matrix(as.numeric(beta), p, 1)
  
  k = p
  k = k - sum(beta==0)
  dummy_beta = beta
  dummy_beta[which(beta==0)] = rnorm(sum(beta==0))
  k = k- sum((dummy_beta[1:(length(dummy_beta)-1)]-dummy_beta[2:length(dummy_beta)])==0)
  
  tau=1
  if(k!=0){
    for(j in 1:k){
      tau = tau*(p-(j-1))/(j);
    }
  }
  gam = 1-log(n)/(2*log(p));
  if(gam<0){print("\n\n\nEBIC_ERROR:gam<0\n\n\n")}
  
  return(-2*sum(dnorm(y,drop(x%*%beta),sqrt(sigma2),log=T))+
           log(n)*k + 2*gam*log(tau))
  
}

##' EBIC calculation.
##' 
##' @param x is a matrix of order n x p where n is number of observations and p is number of predictor variables.
##' @param y y is a vector of response variable of order n x 1.
##' @param beta Regression coefficient.
##' @param sigma2 Dispersion parameter.
##' @param type Estimation method from ex posterior distribution sample.
##' 
##' @importFrom dplyr %>%
##' 
##' @export
ebic = function(x, y, beta, sigma2, type="fusedlasso"){
  if(ncol(x)!=nrow(beta)){
    warning("Incorrect shape of matrix of x and beta : ncol(x)!=nrow(beta) ")
  }
  if(type=="fusedlasso"){
    value=ebic_fl(x, y, beta, sigma2)
  }
  return(value)
}
