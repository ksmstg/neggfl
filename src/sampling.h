#ifndef SAMPLING_H
#define SAMPLING_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double rand_gamma( double theta, double kappa );
double rand_invgauss(double mu, double lambda);

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double genrand_real3(void);

#endif