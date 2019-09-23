rm(list = ls())

library(dplyr)

# 0. Import packages ------------------------------------------------------

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

# 4. Execute --------------------------------------------------------------

mod <- bfl(x, y, lambda2=500)
plot(beta,col="blue",type="b",pch=1,ylim=range(beta, mod$beta))
lines(mod$beta, type="b",lty=1,col="black")
legend("topright",pch=1,lty=1,merge=TRUE,text.col=c("blue","black"),legend=c("True","Fitted"))

lambda2s <- c(50, 100, 150, 250, 500)

hparams <- expand.grid(lambda2=lambda2s)
hparams <- rbind(hparams, 1)
hparams <- rbind(hparams, 100)

grid_res <- list()
grid_ebic <- rep(NA, nrow(hparams))
for(k in 1:nrow(hparams)){
  lambda2 = hparams$lambda2[k]
  mod <- bfl(x, y, lambda2=lambda2)
  
  grid_res[[k]] <- list(beta = mod$beta,  lambda2 = lambda2)
  grid_ebic[k] <- ebic(x, y, mod$beta, mod$sigma2)
  
  fname = paste0("./debug/result/", "bfl_plot_lambda2=", lambda2, ".png")
  png(fname, width = 300, height = 400)  # 描画デバイスを開く
  plot(beta,col="blue",type="b",pch=1,ylim=range(beta, mod$beta))
  lines(mod$beta, type="b",lty=1,col="black")
  legend("topright",pch=1,lty=1,merge=TRUE,text.col=c("blue","black"),legend=c("True","Fitted"))
  dev.off()                      
}

res <- grid_res[[grid_ebic %>% which.max()]]

  
