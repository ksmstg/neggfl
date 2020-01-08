library(Rcpp)

sourceCpp("./src/gbsGflasso.cpp")
sourceCpp("./src/kernelDensity.cpp")
source("./R/sparseAlgorithm.R")

#' @useDynLib neggfl, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

##' Bayesian fused lasso
##' 
##' @param x is a matrix of order n x p where n is number of observations and p is number of predictor variables.
##' @param y y is a vector of response variable of order n x 1.
##' @param lambda2 The hyper-parameter control the degree of fusion.
##' @param maxiter Gibbs Sampling Number of Iterations.
##' @param burnin Number of burnin.
##' @param scaled Flag whether to standardize.
##' @param method Estimation method from ex posterior distribution sample.
##' 
##' @importFrom dplyr %>%
##' @importFrom mvtnorm dmvnorm
##' 
##' @export
bfl <- function(x, y, lambda2, maxiter=1e5, burnin=3e4, scaled=T, method="mode_kde"){
  x_org <- x
  y_org <- y
  
  if(scaled){
    x_sd = x_org %>% var %>% diag %>% sqrt
    x <- scale(x_org)
    y <- y - mean(y_org)
  }
  
  hparams = data.frame(lamb2=lambda2)
  B = maxiter
  posterior = bfl_bas(as.matrix(y), as.matrix(x), hparams, B)
  
  beta_hat_tmp = posterior$beta[,,burnin:B] %>% apply(., 1, eval(parse(text = method))) %>% matrix()
  sigma2_hat = posterior$sigma[,,burnin:B] %>% eval(parse(text = method))()
  
  sfa_res <- sfa_bfl(y, x, beta_hat_tmp, sigma2_hat, lambda2)
  beta_hat <- sfa_res$beta %>% matrix()
  
  if(scaled){
    beta_hat_tmp <-  diag(1/x_sd) %*% beta_hat_tmp
    beta_hat <- diag(1/x_sd) %*% beta_hat
  }
  
  return(list(
    name="bfl",
    beta=beta_hat, 
    beta_index=sfa_res$index, 
    beta_tmp=beta_hat_tmp, 
    sigma2=sigma2_hat,
    lambda2=lambda2))
}
