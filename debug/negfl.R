rm(list = ls())

library(dplyr)

# 0. Import packages ------------------------------------------------------

# library(Rcpp)
# library(RcppEigen)
# library(RcppArmadillo)
# library(BH)
# library(GIGrvg)
# 
# sourceCpp("./src/gbsNegGflasso.cpp")
# library(mvtnorm)
# library(dplyr)
# 
# source("./R/gbsNegGflasso.R")

library("neggfl")

# 1. Data setting ---------------------------------------------------------

n=300
p=100

# 2. Generate data. -------------------------------------------------------

beta=rep(0,p)
beta[1:20]=1
beta[11:15]=2
beta[25]=3
beta[41:45]=1
x=matrix(rnorm(n*p),n,p)
y=x%*%beta+rnorm(n,0,0.5)

# 3. Set hyperparameter. --------------------------------------------------

lambda2=20000.0
gamma2=0.025

# 4. Execute --------------------------------------------------------------

mod <- negfl(x, y, lambda2=3000, gamma2=0.4, maxiter=10000, burnin=3000)
plot(beta,col="blue",type="b",pch=1,ylim=range(beta, mod$beta))
lines(mod$beta, type="b",lty=1,col="black")
legend("topright",pch=1,lty=1,merge=TRUE,text.col=c("blue","black"),legend=c("True","Fitted"))

lambda2s <- c(2500, 2750, 3000, 3250, 3500)
gamma2s <- c(0.15, 0.2, 0.25)

hparams <- expand.grid(lambda2=lambda2s, gamma2=gamma2s)
hparams <- rbind(hparams, c(100, 0.2))
hparams <- rbind(hparams, c(10000, 0.2))

grid_res <- list()
grid_ebic <- rep(NA, nrow(hparams))
for(k in 1:nrow(hparams)){
  hparam = hparams[k,]
  lambda2 = hparam$lambda2
  gamma2 = hparam$gamma2
  mod <- negfl(x, y, lambda2=lambda2, gamma2=gamma2)
  
  grid_res[[k]] <- list(beta = mod$beta,  lambda2 = lambda2, gamma2 = gamma2)
  grid_ebic[k] <- ebic(x, y, mod$beta, mod$sigma2)
  
  fname = paste0("./debug/result/", "plot_lambda2=", lambda2, "_gamma2=", gamma2, ".png")
  png(fname, width = 300, height = 400)  # 描画デバイスを開く
  plot(beta,col="blue",type="b",pch=1,ylim=range(beta, mod$beta))
  lines(mod$beta, type="b",lty=1,col="black")
  legend("topright",pch=1,lty=1,merge=TRUE,text.col=c("blue","black"),legend=c("True","Fitted"))
  dev.off()                      
}

res <- grid_res[[grid_ebic %>% which.max()]]

  
