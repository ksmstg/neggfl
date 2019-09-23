// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "parabolicCylinder.h"

double D(double w, double x, double lamb, double a){
  /* Used report */
  return std::pow(w,(2.0*lamb+a-1.0))*std::exp(-0.5*w*w-x*w);
}

// [[Rcpp::export]]
double Int_simt_log_D(double x, double lamb, double a){
  /* Used report */
  int k, k0, N=100000;
  double s, f, w = 0.00, dh=0.01, DELTA = 1.0e-10;
  s = D(w, x, lamb, a)/2.0;
  
  for (k0=1; k0<=100; k0++){
    w = 0.00;
    dh = 0.001/(k0*k0);
    s = D(w, x, lamb, a)/2.0;
    for (k=1; k<=N; k++) {
      w = w + dh;
      f = D(w, x, lamb, a);
      s = s + f;
      if (std::abs(f) < DELTA) {
        if(dh*s!=0){
          if(std::isinf(std::abs(std::log(dh*s) - x*x/4.0))){
            std::cout << "/n/n/n/n/n logD is Inf /n/n/n/n/n"  << "\n";
          }
          return std::log(dh*s) - x*x/4.0;
        }
      }
    }
  }
  //printf("NOT CONVERGENT\n");
  return DELTA;
}
