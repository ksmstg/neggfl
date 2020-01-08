library(Rcpp)
library(mvtnorm)

sourceCpp("./src/parabolicCylinder.cpp")

dbfl <- function(y, x, beta_hat, sigma2_hat, lamb2){
  score = dmvnorm(as.numeric(y), as.numeric(x%*%beta_hat), sigma2_hat*diag(n), log=T)
  
  dprior_log <- function(beta_hat, sigma2_hat, lamb2){
    p <- length(beta_hat)
    score = 0
    for(j in 2:p){
      score = score - lamb2*abs(beta_hat[j-1]-beta_hat[j])/sqrt(sigma2_hat)
    }
    return(score)
  }
  
  score = score + dprior_log(beta_hat, sigma2_hat, lamb2)
  
  return(score)
}


sfa_bfl <- function(y, x, beta_hat, sigma2_hat, lamb2, repeat_num=5){
  n <- nrow(x)
  p <- ncol(x)
  
  ind_beta = 1:p
  beta_res = beta_hat
  
  scores = c()
  
  for(k in 1:repeat_num){
    for(j in 1:p){
      beta_front = beta_res
      beta_back = beta_res
      
      scores["orig"] = dbfl(y, x, beta_res, sigma2_hat, lamb2)
      
      if(j <= 1){
        scores["front"] = log(0)
      }else if( ind_beta[j-1] != -1 & ind_beta[j-1] != ind_beta[j]){
        beta_front[ind_beta == ind_beta[j]] = beta_res[j-1]
        scores["front"] = dbfl(y, x, beta_front, sigma2_hat, lamb2)
      }else{
        scores["front"] = log(0)
      }
      
      if(j >= p){
        scores["back"] = log(0)
      }else if(ind_beta[j+1] != -1 & ind_beta[j] != ind_beta[j+1]){
        beta_back[ind_beta == ind_beta[j]] = beta_res[j+1]
        scores["back"] = dbfl(y, x, beta_back, sigma2_hat, lamb2)
      }else{
        scores["back"] = log(0)
      }
      
      scores_max = which.max(scores)
      if( names(scores_max) == "front" & j > 1){
        beta_res = beta_front;
        ind_beta[ind_beta == ind_beta[j]] = ind_beta[j-1];
      }else if( names(scores_max) == "back" & j < p){
        beta_res = beta_back;
        ind_beta[ind_beta == ind_beta[j]] = ind_beta[j+1];
      }
    }
  }
  
  return(list(beta=beta_res, index=ind_beta))
  
}

dnegfl <- function(y, x, beta_hat, sigma2_hat, lamb2, gamma2){
  n = nrow(x)
  score = dmvnorm(as.numeric(y), as.numeric(x%*%beta_hat), sigma2_hat*diag(n), log=T)
  
  dprior_log <- function(beta_hat, sigma2_hat, lamb2, gamma2){
    p <- length(beta_hat)
    score = 0
    for(j in 2:p){
      score = score - (beta_hat[j-1]-beta_hat[j])^2 / (4.0*gamma2^2*sigma2_hat) +
        Int_simt_log_D(abs(beta_hat[j-1]-beta_hat[j])/(gamma2*sqrt(sigma2_hat)), lamb2, 1.0)
    }
    return(score)
  }
  
  score = score + dprior_log(beta_hat, sigma2_hat, lamb2, gamma2)
  
  return(score)
}

sfa_negfl <- function(y, x, beta_hat, sigma2_hat, lamb2, gamma2, repeat_num=5){
  n <- nrow(x)
  p <- ncol(x)

  ind_beta = 1:p
  beta_res = beta_hat

  scores = c()

  for(k in 1:repeat_num){
    for(j in 1:p){
      beta_front = beta_res
      beta_back = beta_res

      scores["orig"] = dnegfl(y, x, beta_res, sigma2_hat, lamb2, gamma2)

      if(j <= 1){
        scores["front"] = log(0)
      }else if( ind_beta[j-1] != -1 & ind_beta[j-1] != ind_beta[j]){
        beta_front[ind_beta == ind_beta[j]] = beta_res[j-1]
        scores["front"] = dnegfl(y, x, beta_front, sigma2_hat, lamb2, gamma2)
      }else{
        scores["front"] = log(0)
      }

      if(j >= p){
        scores["back"] = log(0)
      }else if(ind_beta[j+1] != -1 & ind_beta[j] != ind_beta[j+1]){
        beta_back[ind_beta == ind_beta[j]] = beta_res[j+1]
        scores["back"] = dnegfl(y, x, beta_back, sigma2_hat, lamb2, gamma2)
      }else{
        scores["back"] = log(0)
      }

      scores_max = which.max(scores)
      if( names(scores_max) == "front" & j > 1){
        beta_res = beta_front;
        ind_beta[ind_beta == ind_beta[j]] = ind_beta[j-1];
      }else if( names(scores_max) == "back" & j < p){
        beta_res = beta_back;
        ind_beta[ind_beta == ind_beta[j]] = ind_beta[j+1];
      }

      # print(c(beta_res))
      # print(ind_beta)
    }
  }

  return(list(beta=beta_res, index=ind_beta))

}

# dnegsfl <- function(y, x, beta_hat, sigma2_hat, lamb2, gamma2){
#   score = dmvnorm(as.numeric(y), as.numeric(x%*%beta_hat), sigma2_hat*diag(n), log=T)
#   
#   dprior_log <- function(y, x, beta_hat, sigma2_hat, lamb2, gamma2){
#     p <- ncol(x)
#     score = 0
#     for(j in 2:p){
#       score = score + (beta_hat[j-1]-beta_hat[j])^2 / (4.0*gamma2^2*sigma2_hat) + 
#         Int_simt_log_D(abs(beta_hat[j-1]-beta_hat[j])/(gamma2*sqrt(sigma2_hat)), lamb2, 1.0)
#     } 
#     return(score)
#   }
#   
#   score = score + dprior_log(y, x, beta_hat, sigma2_hat, lamb2, gamma2)
#   
#   return(score)
# }

# sfa_negsfl <- function(y, x, beta_hat, sigma2_hat, lamb2, gamma2, repeat_num=5){
#   n <- nrow(x)
#   p <- ncol(x)
#   
#   ind_beta = 1:p
#   beta_res = beta_hat
#   
#   scores = c()
#   
#   for(k in 1:repeat_num){
#     for(j in 1:p){
#       beta_sparse = beta_res
#       beta_front = beta_res
#       beta_back = beta_res
#       
#       scores["orig"] = dnegsfl(y, x, beta_res, sigma2_hat, lamb2, gamma2)
#       
#       beta_sparse[ind_beta == ind_beta[j]] = 0
#       scores["sparse"] = dnegfl(y, x, beta_sparse, sigma2_hat, lamb2, gamma2)
#       
#       if(j <= 1){
#         scores["front"] = log(0)
#       }else if( ind_beta[j-1] != -1 & ind_beta[j-1] != ind_beta[j]){
#         beta_front[ind_beta == ind_beta[j]] = beta_res[j-1]
#         scores["front"] = dnegfl(y, x, beta_front, sigma2_hat, lamb2, gamma2)
#       }else{
#         scores["front"] = log(0)
#       }
#       
#       if(j >= p){
#         scores["back"] = log(0)
#       }else if(ind_beta[j+1] != -1 & ind_beta[j] != ind_beta[j+1]){
#         beta_back[ind_beta == ind_beta[j]] = beta_res[j+1]
#         scores["back"] = dnegfl(y, x, beta_back, sigma2_hat, lamb2, gamma2)
#       }else{
#         scores["back"] = log(0)
#       }
#       
#       scores_max = which.max(scores)
#       if( names(scores_max) == "sparse" ){
#         beta_res[ind_beta == ind_beta[j]] = 0.0;
#         ind_beta[ind_beta == ind_beta[j]] = -1;
#       }else if( names(scores_max) == "front" & j > 1){
#         beta_res[ind_beta == ind_beta[j]] = beta_res[j-1];
#         ind_beta[ind_beta == ind_beta[j]] = ind_beta[j-1];
#       }else if( names(scores_max) == "back" & j < p){
#         beta_res[ind_beta == ind_beta[j]] = beta_res[j+1];
#         ind_beta[ind_beta == ind_beta[j]] = ind_beta[j+1];
#       }
#     }
#   }
#   
#   return(beta_res)
#   
# }
