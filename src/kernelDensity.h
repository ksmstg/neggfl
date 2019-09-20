#ifndef KERNELDENSITY_H
#define KERNELDENSITY_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double bw_nrd0(arma::vec x);
float Ko(double t);
double mode_kde(arma::vec points, int T);

#endif