// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "utils.h"

double sq(double x){
  return x*x;
}

double sq_norm(arma::vec x){
  double tmp = norm(x, 2);
  return tmp*tmp;
}

double sq_norm_mat(arma::mat x){
  double tmp = norm(x, 2);
  return tmp*tmp;
}



