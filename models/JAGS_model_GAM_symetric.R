## model used to smooth series of annual indices, that have an estimated precision

model {
  
  for (i in 1:ncounts) { 
    lindex[i] ~ dnorm(mu[i],preci[i]) # response
    mu[i] ~ dnorm(x.gam[i],tau) #expected response
    
  } ## expected response
  
  sd <- 1/pow(tau,0.5)
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/10^2 is appropriate!
  
  
  x.gam <- X %*% b ## expected response
  
  for (i in 1:1) { b[i] ~ dnorm(0,0.01) }
  ## prior for s(year)... 
  K1 <- S1[1:(nknots-1),1:(nknots-1)] * lambda[1]  + S1[1:(nknots-1),(nknots):((nknots*2)-2)] * lambda[2]
  b[2:(nknots)] ~ dmnorm(zero[2:(nknots)],K1) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  
  
  
  
  
  ### derived parameters
  ind.pred.l <- X.pred %*% b #full predicted values
  for(y in 1:nyears){
    ind.pred[y] <- exp(ind.pred.l[y])
  }
  Ysb <- max(Ys,Yb-(Ye-Yb))#start point for change in trend comparison
  C1 <- ind.pred[Yb]/ind.pred[Ysb] #proportional change from Ysb to Yb - early short-term
  C2 <- ind.pred[Ye]/ind.pred[Yb] #proportional change from Yb to Ye - recent short-term
  C3 <- ind.pred[Ye]/ind.pred[Ys] #proportional change from Ys to Ye - Long-term change
  T1 <- ((C1)^(1/(Ye-Yb))-1)*100 ##average annual percent change from Ysb to Yb
  T2 <- ((C2)^(1/(Ye-Yb))-1)*100 ##average annual percent change from Yb to Ye
  T3 <- ((C3)^(1/(Ye-Ys))-1)*100 ##average annual percent change from Ys to Ye - long term change
  
  Tdif <- T2-T1
  Tdif_neg = step(-1*Tdif)
  
  
  
}