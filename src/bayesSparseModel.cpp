// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include "bayesSparseModel.h"

bsm::bsm(arma::mat _y, arma::mat _X, int _B):
  n(_X.n_rows), p(_X.n_cols), B(_B), y(_y), X(_X), In(arma::eye(n,n)) {
  XtX = trans(X)*X;
  Xy = trans(X)*y;
  
  nu0_sigma2 = 0.00001, eta0_sigma2 = 0.00001;
  
  beta_list = arma::cube(p,1,B).fill(NA_REAL);
  beta_list.slice(0) =  arma::zeros(p,1);
  
  sigma2_list = arma::cube(1,1,B).fill(NA_REAL);
  sigma2_list.slice(0) = 1.0;
}

bsm::~bsm(){};

void bsm::up(Rcpp::List params, int t){
  std::cout << "message from bsm" << std::endl;
}

Rcpp::List bsm::gibbs_sampler() {
  if(hparams.nrows()==1){
    Progress prober(B-1, true);
    for(t=0;t<(B-1);t++){
      prober.increment(); 
      up(params, t);
    }
    return params;
  } else{
    // TODO:複数の超パラメータに対する処理を追加する
    return params;
  }
}
