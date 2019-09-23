#ifndef GBSGFLASSO_H
#define GBSGFLASSO_H

// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "bayesSparseModel.h"

// Bayes fuzed lasso modeling
class bflm : public bsm{
protected:
  arma::mat invtau2_til;
  arma::mat _invtau2_til;
  arma::mat invtau2_til_new;
  arma::mat invtau2_til_old;
  double lamb_invtau2_til, mu_invtau2_til;
  arma::cube invtau2_til_list;
  
  double lamb2;
  
public:
  bflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B);
  ~bflm();
};

class bflm_bas : public bflm{
public:
  bflm_bas(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B);
  ~bflm_bas();
  
  arma::mat make_S_beta(arma::mat _invtau2_til);
  
  arma::mat samp_beta(arma::mat _invtau2_til, double _sigma2);
  arma::mat samp_invtau2_til(arma::mat _beta, double _sigma2);
  double samp_sigma2(arma::mat _beta, arma::mat _invtau2_til);
  
  void up(Rcpp::List params, int t);
};

Rcpp::List bfl_bas(arma::mat y, arma::mat X, Rcpp::DataFrame hparams, int B);

#endif