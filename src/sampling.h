#ifndef SAMPLING_H
#define SAMPLING_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double rand_gamma( double theta, double kappa );
double rand_invgauss(double mu, double lambda);

#endif