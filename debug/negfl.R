rm(list = ls())

# 0. Import packages ------------------------------------------------------

library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
library(BH)
library(GIGrvg)

sourceCpp("./src/gbsNegGflasso.cpp")
library(mvtnorm)
library(dplyr)

source("./R/gbsNegGflasso.R")


# 1. Data setting ---------------------------------------------------------

n=50
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

mod <- negfl(x, y, lambda2=20000.0, gamma2=0.025, maxiter=10000, burnin=3000)
plot(beta,col="blue",type="b",pch=1,ylim=range(beta, mod$beta))
lines(mod$beta, type="b",lty=1,col="black")
legend("topright",pch=1,lty=1,merge=TRUE,text.col=c("blue","black"),legend=c("True","Fitted"))

