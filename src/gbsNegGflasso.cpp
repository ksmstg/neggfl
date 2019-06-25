#include "sampling.h"
// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppDist.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


double sq(double x){
  return x*x;
}

double sq_norm(arma::vec x){
  double tmp = norm(x, 2);
  return tmp*tmp;
}

double sq_norm_mat(arma::mat x){
  double tmp = norm(x, 2);
  return tmp*tmp;
}

double sqrt(double);
  
double rinvgauss(double mu, double lambda){
  
  double z,y,x,u, random;
  
  z=R::rnorm(0,1);
  y=z*z;
  x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
  u=R::runif(0,1);
  if(u <= mu/(mu+x)){
    random = x;
  }else{
    random = mu*mu/x;
  };
  
  return(random);
}

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
    //dh = (0.001+std::sqrt(std::abs(f)))*k0;
    for (k=1; k<=N; k++) {
      w = w + dh;
      f = D(w, x, lamb, a);
      //if(k==1){std::cout << k << "\t" << "f : "<< f << "\n";}
      s = s + f;
      //std::cout << "test01 : "<< k0 << "\t" << k << "\t" << std::abs(f) << "\n";
      if (std::abs(f) < DELTA) {
        if(dh*s!=0){
          //std::cout << "CONVERGENT"<< dh*s*std::exp(-x*x/4.0) << "\n";
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

// Bayes sparse modeling class
class bflm{
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
    bflm(arma::mat _y, arma::mat _X, int _B):
    n(_X.n_rows), p(_X.n_cols), B(_B), y(_y), X(_X), In(arma::eye(n,n)) {
      XtX = trans(X)*X;
      Xy = trans(X)*y;
        
      nu0_sigma2 = 0.00001, eta0_sigma2 = 0.00001;
      
      beta_list = arma::cube(p,1,B).fill(NA_REAL);
      beta_list.slice(0) =  arma::zeros(p,1);

      sigma2_list = arma::cube(1,1,B).fill(NA_REAL);
      sigma2_list.slice(0) = 1.0;
    }
    
    virtual ~bflm() {}
    
    virtual void up(Rcpp::List params, int t){
      std::cout << "message from bsm" << std::endl;
    }
    
    Rcpp::List gibbs_sampler() {
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
};

// Bayes Neg fuzed lasso modeling
class negflm : public bflm{
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
  negflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B): 
  bflm(_y, _X, _B){
    
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
  }
  
  ~negflm() {}
};

class negflm_bas : public negflm{
  public:
    negflm_bas(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B) : 
    negflm(_y, _X, _hparams, _B){}
    
    ~negflm_bas() {}
    
    arma::mat make_S_beta(arma::mat _invtau2_til){
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
    }
    
    arma::mat samp_beta(arma::mat _invtau2_til, double _sigma2){
      S_beta = make_S_beta(_invtau2_til);
      
      invS_beta = inv(XtX+S_beta);
      beta_new = rmvnorm(1, invS_beta*Xy, _sigma2*invS_beta).t();

      return beta_new;
    }
    
    arma::mat samp_invtau2_til(arma::mat _beta, arma::mat _psi_til, double _sigma2){
      invtau2_til_new = arma::zeros(p-1,1);
      for(int j=0 ; j < p-1 ; j++){
        lamb_invtau2_til = 2.0 * _psi_til(j,0);
        mu_invtau2_til = sqrt( (lamb_invtau2_til*_sigma2) / std::abs(_beta(j+1,0)-_beta(j,0)) );
        invtau2_til_new(j,0) = rinvgauss(mu_invtau2_til, lamb_invtau2_til);
      }
      return invtau2_til_new;
    }
    
    arma::mat samp_psi_til(arma::mat _invtau2_til, double _sigma2){
      psi_til_new = arma::zeros(p-1,1);
      for(int j=0 ; j < p-1 ; j++){
        theta_psi_til = _invtau2_til(j,0) / (1.0 + _invtau2_til(j,0)*sq(gamma2));
        kappa_psi_til = lamb2 + 1.0;
        psi_til_new(j,0) = rand_gamma(theta_psi_til, kappa_psi_til);
      }
      return psi_til_new;
    }
    
    double samp_sigma2(arma::mat _beta, arma::mat _invtau2_til){
      S_beta = make_S_beta(_invtau2_til);

      eta1_sigma2 = ( sq_norm_mat(y-X*_beta) + accu(trans(_beta)*S_beta*_beta) + eta0_sigma2 )/2.0;
      nu1_sigma2 = ( double(n+p-1) + nu0_sigma2 )/2.0;
      sigma2_new = 1.0 / rand_gamma(1/eta1_sigma2,nu1_sigma2);
      return  sigma2_new;
    }
  
    // double post_weight(Rcpp::List param){
    //   
    //   beta        = arma::as<mat>(param["beta"]);
    //   invtau2_til = arma::as<mat>(param["invtau2_til"]);
    //   psi_til     = arma::as<mat>(param["psi_til"]);
    //   sigma2      = arma::as<double>(param["sigma2"]);
    //   
    //   lamb2 = hparams["lamb2"];
    //   gamma2 = hparams["gamma2"];
    //   
    //   double score_X = (-double(n*p)/2.0)*log(2.0*M_PI*sigma2) + ( -sqNorm_mat(X-A)/(2.0*sigma2) );
    //   
    //   // double score_A = (-double(p*w_len)/2.0) * log(2.0*M_PI*sigma2);
    //   double score_A = -0.5 * double(w_len) * log(sigma2);
    //   double score_invtau2 = 0.0;
    //   for(int k=0 ; k < w_len ; k++){
    //     i1 = w_ind(k,0);
    //     i2 = w_ind(k,1);
    //     
    //     score_A       += 0.5 * log( invtau2(i1,i2) ) + ( -sqNorm_mat(A.row(i1)-A.row(i2)) / (2.0*sigma2) * invtau2(i1,i2) );
    //     score_invtau2 += log(lamb*lamb*w(i1,i2)*w(i1,i2)/2.0) + ( -lamb*lamb*w(i1,i2)*w(i1,i2)/2.0/invtau2(i1,i2) );
    //   }
    //   
    //   double score_sigma2 = (-nu0_sigma2/2.0+1) * log(sigma2) + (-eta0_sigma2/2.0/sigma2);
    //   
    //   return score_X + score_A + score_invtau2 + score_sigma2 ;
    // }
    
    void up(Rcpp::List params, int t) override {
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
    }
};


// Bayes Neg sparse fuzed lasso modeling
class negsflm : public bflm{
protected:
  arma::mat invtau2_new;
  arma::mat invtau2_old;
  double lamb_invtau2, mu_invtau2;
  arma::mat invtau2_til_new;
  arma::mat invtau2_til_old;
  double lamb_invtau2_til, mu_invtau2_til;
  
public:
  negsflm(arma::mat _y, arma::mat _X, Rcpp::DataFrame _hparams, int _B): 
  bflm(_y, _X, _B){
    
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
  }
  
  ~negsflm() {}
};


// [[Rcpp::export]]
Rcpp::List negfl_bas(arma::mat y, arma::mat X, Rcpp::DataFrame hparams, int B){
  negflm_bas model(y, X, hparams, B);
  return model.gibbs_sampler();
}

