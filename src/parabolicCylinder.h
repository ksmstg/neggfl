#ifndef PARABOLICCYLINDER_H
#define PARABOLICCYLINDER_H

// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

double D(double w, double x, double lamb, double a);
double Int_simt_log_D(double x, double lamb, double a);

#endif