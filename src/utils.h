#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double sq(double x);
double sq_norm(arma::vec x);
double sq_norm_mat(arma::mat x);
double sqrt(double);

#endif