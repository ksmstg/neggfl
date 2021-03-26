// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

#include "gbsNegGflasso.h"
#include "bayesSparseModel.h"
#include "sampling.h"
#include "utils.h"

negflm::negflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B): 
  bsm(_y, _X, _B){
    
  invtau2_til_list = arma::cube(p-1,1,B).fill(NA_REAL);
  invtau2_til_list.slice(0) = arma::ones(p-1,1);
  psi_til_list = arma::cube(p-1,1,B).fill(NA_REAL);
  psi_til_list.slice(0) = arma::ones(p-1,1);
  
  hparams = _hparams;
  gamma2 = hparams["gamma2"];
  lamb2 = hparams["lamb2"];
  
  params = Rcpp::List::create( 
    Rcpp::Named("beta") = beta_list,
    Rcpp::Named("invtau2_til") = invtau2_til_list,
    Rcpp::Named("psi_til") = psi_til_list,
    Rcpp::Named("sigma2")  = sigma2_list
  );
};
  
negflm::~negflm() {};

negflm_bas::negflm_bas(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B) : 
  negflm(_y, _X, _hparams, _B){};
    
negflm_bas::~negflm_bas() {};

arma::mat negflm_bas::make_S_beta(arma::mat _invtau2_til){
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
  
arma::mat negflm_bas::samp_beta(arma::mat _invtau2_til, double _sigma2){
  S_beta = make_S_beta(_invtau2_til);
  
  invS_beta = inv(XtX+S_beta);
  invS_beta = (invS_beta+invS_beta.t())/2;
  beta_new = rmvnorm(1, invS_beta*Xy, _sigma2*invS_beta).t();

  return beta_new;
};
  
arma::mat negflm_bas::samp_invtau2_til(arma::mat _beta, arma::mat _psi_til, double _sigma2){
  invtau2_til_new = arma::zeros(p-1,1);
  for(int j=0 ; j < p-1 ; j++){
    lamb_invtau2_til = 2.0 * _psi_til(j,0);
    mu_invtau2_til = sqrt( (lamb_invtau2_til*_sigma2) ) / std::abs(_beta(j+1,0)-_beta(j,0));
    invtau2_til_new(j,0) = rand_invgauss(mu_invtau2_til, lamb_invtau2_til);
  }
  return invtau2_til_new;
};
    
arma::mat negflm_bas::samp_psi_til(arma::mat _invtau2_til, double _sigma2){
  psi_til_new = arma::zeros(p-1,1);
  for(int j=0 ; j < p-1 ; j++){
    theta_psi_til = _invtau2_til(j,0) / (1.0 + _invtau2_til(j,0)*sq(gamma2));
    kappa_psi_til = lamb2 + 1.0;
    psi_til_new(j,0) = rand_gamma(theta_psi_til, kappa_psi_til);
  }
  return psi_til_new;
};
    
double negflm_bas::samp_sigma2(arma::mat _beta, arma::mat _invtau2_til){
  S_beta = make_S_beta(_invtau2_til);

  eta1_sigma2 = ( sq_norm_mat(y-X*_beta) + accu(trans(_beta)*S_beta*_beta) + eta0_sigma2 )/2.0;
  nu1_sigma2 = ( double(n+p-1) + nu0_sigma2 )/2.0;
  sigma2_new = 1.0 / rand_gamma(1/eta1_sigma2,nu1_sigma2);
  return  sigma2_new;
};

void negflm_bas::up(Rcpp::List params, int t) {
  _beta        = Rcpp::as<arma::cube>(params["beta"]).slice(t);
  _invtau2_til = Rcpp::as<arma::cube>(params["invtau2_til"]).slice(t);
  _psi_til     = Rcpp::as<arma::cube>(params["psi_til"]).slice(t);
  _sigma2      = Rcpp::as<arma::cube>(params["sigma2"]).slice(t)(0,0);
  
  _beta        = samp_beta(_invtau2_til, _sigma2);
  _invtau2_til = samp_invtau2_til(_beta, _psi_til, _sigma2);
  _psi_til     = samp_psi_til(_invtau2_til, _sigma2);
  _sigma2      = samp_sigma2(_beta, _invtau2_til);

  // logLike_new   = culc_logLike(A_new, sigma2_new);
  
  Rcpp::as<arma::cube>(params["beta"]).slice(t+1)        = _beta;
  Rcpp::as<arma::cube>(params["invtau2_til"]).slice(t+1) = _invtau2_til;
  Rcpp::as<arma::cube>(params["psi_til"]).slice(t+1)     = _psi_til;
  Rcpp::as<arma::cube>(params["sigma2"]).slice(t+1)      = _sigma2;
  
  // param = List::create( 
  //   Named("A") = A_new,
  //   Named("invtau2") = invtau2_new,
  //   Named("sigma2")  = sigma2_new
  // );
  
  // as<cube>(params["pw_list"]).slice(t+1) = post_weight( A_new, sigma2_new, invtau2_new, hparams);
  // as<cube>(params["pw_list"]).slice(t+1) = post_weight(param);
  // as<cube>(params["logLike_list"]).slice(t+1) = logLike_new;
  
  // if( NumericVector::is_na(sigma2_new) ){
  //   stop("Error: Unexpected condition occurred");
  // }
  // 
  // std::cout << "A_new : " << A_new << std::endl;
  // std::cout << "invtau2_new : " << invtau2_new << std::endl;
  // std::cout << "sigma2_new : " << sigma2_new << std::endl;
};


negsflm::negsflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B): 
  bsm(_y, _X, _B){
  
  arma::cube invtau2_list = arma::cube(p,1,B).fill(NA_REAL);
  invtau2_list.slice(0) = arma::mat(p,1).fill(1.0/p);
  
  arma::cube invtau2_til_list = arma::cube(p,1,B).fill(NA_REAL);
  invtau2_til_list.slice(0) = arma::mat(p,1).fill(1.0/p);
  
  params = Rcpp::List::create( 
    Rcpp::Named("beta_list") = beta_list,
    Rcpp::Named("invtau2_list") = invtau2_list,
    Rcpp::Named("invtau2_til_list") = invtau2_til_list,
    Rcpp::Named("sigma2_list")  = sigma2_list
  );
};
  
negsflm::~negsflm() {};


// [[Rcpp::export]]
Rcpp::List negfl_bas(arma::mat y, arma::mat X, Rcpp::DataFrame hparams, int B){
  negflm_bas model(y, X, hparams, B);
  return model.gibbs_sampler();
};

