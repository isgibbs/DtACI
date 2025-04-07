vecZeroMax <- Vectorize(function(x){max(x,0)})
vecZeroMin <- Vectorize(function(x){min(x,0)})

### Definition of the pinball loss function
pinball <- function(u,alpha){
  alpha*u - vecZeroMin(u)
}

### AgACI method
agaci <- function(betas,alpha,gammas,alphaInit=alpha,eps=0.001){
  T <- length(betas)
  k <- length(gammas)
  alphaSeq <- rep(alphaInit,T)
  errSeq <- rep(0,T)
  gammaSeq <- rep(0,T)

  expertAlphas <- rep(alphaInit,k)
  expertProbs <- rep(1/k,k)
  expertSqLosses <- rep(0,k)
  expertEtas <- rep(0,k)
  expertLValues <- rep(0,k)
  expertMaxLosses <- rep(0,k)
  for(t in 1:T){
    ### compute predictions
    alphaSeq[t] <- sum(expertProbs*expertAlphas)
    errSeq[t] <- as.numeric(alphaSeq[t]>betas[t])
    gammaSeq[t] <- sum(expertProbs*gammas)
    
    ### update expert weights
    expertLosses <- (errSeq[t]-alpha)*(expertAlphas - alphaSeq[t])
    expertSqLosses <- expertSqLosses + expertLosses^2
    expertMaxLosses <- pmax(expertMaxLosses,abs(expertLosses))
    expertEVals <- 2^(ceiling(log(abs(expertMaxLosses)+eps,2)) + 1)
    expertLValues <- expertLValues + (1/2)*(expertLosses*(1+expertEtas*expertLosses) + expertEVals*(expertEtas*expertLosses > 1/2))
    expertEtas <- pmin(1/expertEVals,sqrt(log(k)/expertSqLosses))
    expertAlphas <- expertAlphas + gammas*(alpha-as.numeric(expertAlphas>betas[t]))
    expertWeights <- expertEtas*exp(-expertEtas*expertLValues + max(expertEtas*expertLValues))
    expertProbs <- expertWeights/sum(expertWeights)

  }
  
  return(list(alphaSeq,errSeq,gammaSeq))
}



