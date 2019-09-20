// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "kernelDensity.h"

double bw_nrd0(arma::vec x){
  /* Used report */
  int i,n;
  double lo;
  arma::vec x_sort;
  
  lo=0.0;
  n = x.size();
  for(i=0;i<n;i++){
    lo += (x(i)-arma::mean(x))*(x(i)-arma::mean(x));
  }
  lo = lo/double(n-1);
  
  return 0.9 * lo * std::pow(double(n),-0.2);
}

float Ko(double t){
  /* Used report */
  return std::exp(-t*t)/(2.0*M_PI);
}

// [[Rcpp::export]]
double mode_kde(arma::vec points, int T=512){
  /* Used report */
  int i,j,K = points.size();
  double result = 0;
  float bw = bw_nrd0(points);
  //float bw = 1.0;
  arma::vec kde_points(T);
  arma::vec kde_densvalue;
  
  
  double bounds = std::abs( arma::max(points)-arma::min(points) );
  //std::cout << bounds << "\n";
  kde_densvalue = arma::zeros(T);
  
  for(i=0; i<T ; i++){
    kde_points(i) =	bounds*( (double)i/(double)(T-1) ) + arma::min(points);
    //std::cout << kde_points(i) << "\n";
    for(j=0; j<K; j++){
      kde_densvalue(i) += Ko((kde_points(i)-points(j))/bw);
    }
  }
  
  for(i=0; i<T; i++){
    if(kde_densvalue(i)==arma::max(kde_densvalue)){
      result = kde_points(i);
    }
  }
  
  return result;
}
