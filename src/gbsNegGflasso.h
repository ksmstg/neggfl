#ifndef GBSNEGGFLASSO_H
#define GBSNEGGFLASSO_H

// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "bayesSparseModel.h"

// Bayes Neg fuzed lasso modeling
class negflm : public bsm{
protected:
  arma::mat invtau2_til;
  arma::mat _invtau2_til;
  arma::mat invtau2_til_new;
  arma::mat invtau2_til_old;
  double lamb_invtau2_til, mu_invtau2_til;
  arma::cube invtau2_til_list;
  
  arma::mat psi_til;
  arma::mat _psi_til;
  arma::mat psi_til_new;
  arma::mat psi_til_old;
  double theta_psi_til, kappa_psi_til;
  arma::cube psi_til_list;
  
  double lamb2, gamma2;
  
public:
  negflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B);
  ~negflm();
};

class negflm_bas : public negflm{
public:
  negflm_bas(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B);
  ~negflm_bas();
  
  arma::mat make_S_beta(arma::mat _invtau2_til);
  
  arma::mat samp_beta(arma::mat _invtau2_til, double _sigma2);
  arma::mat samp_invtau2_til(arma::mat _beta, arma::mat _psi_til, double _sigma2);
  arma::mat samp_psi_til(arma::mat _invtau2_til, double _sigma2);
  double samp_sigma2(arma::mat _beta, arma::mat _invtau2_til);
  
  void up(Rcpp::List params, int t);
};


// Bayes Neg sparse fuzed lasso modeling
class negsflm : public bsm{
protected:
  arma::mat invtau2_new;
  arma::mat invtau2_old;
  double lamb_invtau2, mu_invtau2;
  arma::mat invtau2_til_new;
  arma::mat invtau2_til_old;
  double lamb_invtau2_til, mu_invtau2_til;
  
public:
  negsflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B);
  ~negsflm();
};

Rcpp::List negfl_bas(arma::mat y, arma::mat X, Rcpp::DataFrame hparams, int B);

#endif