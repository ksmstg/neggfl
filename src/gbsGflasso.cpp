// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "gbsGflasso.h"
#include "bayesSparseModel.h"
#include "sampling.h"
#include "utils.h"

bflm::bflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B): 
  bsm(_y, _X, _B){
    
  invtau2_til_list = arma::cube(p-1,1,B).fill(NA_REAL);
  invtau2_til_list.slice(0) = arma::ones(p-1,1);

  hparams = _hparams;
  lamb2 = hparams["lamb2"];
  
  params = Rcpp::List::create( 
    Rcpp::Named("beta") = beta_list,
    Rcpp::Named("invtau2_til") = invtau2_til_list,
    Rcpp::Named("sigma2")  = sigma2_list
  );
};
  
bflm::~bflm() {};

bflm_bas::bflm_bas(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B) : 
  bflm(_y, _X, _hparams, _B){};
    
bflm_bas::~bflm_bas() {};

arma::mat bflm_bas::make_S_beta(arma::mat _invtau2_til){
  S_beta = arma::zeros(p,p);
  
  S_beta(0,0) = _invtau2_til(0,0);
  S_beta(0,1) = -1.0 * _invtau2_til(0,0);
  for(int j=1 ; j<p-1 ; j++){
    S_beta(j,j-1) = -1.0 * _invtau2_til(j-1,0);
    S_beta(j,j) = _invtau2_til(j-1,0) + _invtau2_til(j,0);
    S_beta(j,j+1) = -1.0 * _invtau2_til(j,0);
  }
  S_beta(p-1,p-2) = -1.0 * _invtau2_til(p-2,0);
  S_beta(p-1,p-1) = _invtau2_til(p-2,0);

  return S_beta;
};
  
arma::mat bflm_bas::samp_beta(arma::mat _invtau2_til, double _sigma2){
  S_beta = make_S_beta(_invtau2_til);
  
  invS_beta = inv(XtX+S_beta);
  beta_new = rmvnorm(1, invS_beta*Xy, _sigma2*invS_beta).t();

  return beta_new;
};
  
arma::mat bflm_bas::samp_invtau2_til(arma::mat _beta, double _sigma2){
  invtau2_til_new = arma::zeros(p-1,1);
  for(int j=0 ; j < p-1 ; j++){
    lamb_invtau2_til = lamb2*lamb2;
    mu_invtau2_til = sqrt( (lamb_invtau2_til*_sigma2) ) / std::abs(_beta(j+1,0)-_beta(j,0));
    invtau2_til_new(j,0) = rand_invgauss(mu_invtau2_til, lamb_invtau2_til);
  }
  return invtau2_til_new;
};

double bflm_bas::samp_sigma2(arma::mat _beta, arma::mat _invtau2_til){
  S_beta = make_S_beta(_invtau2_til);

  eta1_sigma2 = ( sq_norm_mat(y-X*_beta) + accu(trans(_beta)*S_beta*_beta) + eta0_sigma2 )/2.0;
  nu1_sigma2 = ( double(n+p-1) + nu0_sigma2 )/2.0;
  sigma2_new = 1.0 / rand_gamma(1/eta1_sigma2,nu1_sigma2);
  return  sigma2_new;
};

void bflm_bas::up(Rcpp::List params, int t) {
  _beta        = Rcpp::as<arma::cube>(params["beta"]).slice(t);
  _invtau2_til = Rcpp::as<arma::cube>(params["invtau2_til"]).slice(t);
  _sigma2      = Rcpp::as<arma::cube>(params["sigma2"]).slice(t)(0,0);
  
  _beta        = samp_beta(_invtau2_til, _sigma2);
  _invtau2_til = samp_invtau2_til(_beta, _sigma2);
  _sigma2      = samp_sigma2(_beta, _invtau2_til);

  // logLike_new   = culc_logLike(A_new, sigma2_new);
  
  Rcpp::as<arma::cube>(params["beta"]).slice(t+1)        = _beta;
  Rcpp::as<arma::cube>(params["invtau2_til"]).slice(t+1) = _invtau2_til;
  Rcpp::as<arma::cube>(params["sigma2"]).slice(t+1)      = _sigma2;
                                                                                                                                                                                                                                                                                                         };

// [[Rcpp::export]]
Rcpp::List bfl_bas(arma::mat y, arma::mat X, Rcpp::DataFrame hparams, int B){
  bflm_bas model(y, X, hparams, B);
  return model.gibbs_sampler();
};

