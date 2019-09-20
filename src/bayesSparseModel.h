#ifndef BAYESSPARSEMODEL_H
#define BAYESSPARSEMODEL_H

// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

// Bayes sparse modeling class
class bsm
{
protected:
  int n, p, t, B;
  arma::mat y;
  arma::mat X;
  arma::mat XtX;
  arma::mat Xy;
  
  arma::mat beta;
  arma::mat _beta;
  arma::mat beta_old;
  arma::mat beta_new;
  arma::cube beta_list;
  arma::mat S_beta, invS_beta;
  
  double sigma2;
  double _sigma2;
  double sigma2_old;
  double sigma2_new;
  double nu0_sigma2, eta0_sigma2, nu1_sigma2, eta1_sigma2;
  arma::cube sigma2_list;
  
  Rcpp::List param;
  Rcpp::List params;
  Rcpp::DataFrame hparams;
  
  arma::mat In;
  
public:
  bsm(arma::mat _y, arma::mat _X, int _B);
  
  virtual ~bsm();
  
  virtual void up(Rcpp::List params, int t);
  
  Rcpp::List gibbs_sampler();
};

#endif