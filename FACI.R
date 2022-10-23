vecZeroMax <- Vectorize(function(x){max(x,0)})
vecZeroMin <- Vectorize(function(x){min(x,0)})

### Definition of the pinball loss function
pinball <- function(u,alpha){
  alpha*u - vecZeroMin(u)
}

### Input values are the sequence beta_t, the target level, and the sequence of candidate gammas.
### Return value is a list containing the vectors alpha_t, err_t(alpha_t), err_t(alpha), 
### gamma_t, alphabar_t, err_t(alphabar_t), gammabar_t. Here we use the notation err_t(x)
### to refer to the errors computed using input x.
conformalAdaptStable <- function(betas,alpha,gammas,sigma=1/500,eta=2.8){
  T <- length(betas)
  k <- length(gammas)
  alphaSeq <- rep(alpha,T)
  errSeqAdapt <- rep(0,T)
  errSeqFixed <- rep(0,T)
  gammaSeq <- rep(0,T)
  meanAlphaSeq <- rep(0,T)
  meanErrSeq <- rep(0,T)
  meanGammas <- rep(0,T)
  
  expertAlphas <- rep(alpha,k)
  expertWs <- rep(1,k)
  curExpert <- sample(1:k,1)
  expertCumulativeLosses <- rep(0,k)
  expertProbs <- rep(1/k,k)
  for(t in 1:T){
    alphat <- expertAlphas[curExpert]
    alphaSeq[t] <- alphat
    errSeqAdapt[t] <- as.numeric(alphat>betas[t])
    errSeqFixed[t] <- as.numeric(alpha>betas[t])
    gammaSeq[t] <- gammas[curExpert]
    meanAlphaSeq[t] <- sum(expertProbs*expertAlphas)
    meanErrSeq[t] <- as.numeric(meanAlphaSeq[t]>betas[t])
    meanGammas[t] <- sum(expertProbs*gammas)
    
    #loss <- pinball(betas[t] - alphat,alpha)
    expertLosses <- pinball(betas[t] - expertAlphas,alpha)
    
    expertAlphas <- expertAlphas + gammas*(alpha-as.numeric(expertAlphas>betas[t]))
    
    if(eta < Inf){
      expertBarWs <- expertWs*exp(-eta*expertLosses)
      #expertNextWs <- (1-sigma)*expertBarWs + sigma*(sum(expertBarWs))/k
      expertNextWs <- (1-sigma)*expertBarWs/sum(expertBarWs) + sigma/k
      
      #switch <- rbinom(1,1,max(1-expertNextWs[curExpert]/expertWs[curExpert],0))
      expertProbs <- expertNextWs/sum(expertNextWs)
      #if(switch == 1){
      curExpert <- sample(1:k,1,prob=expertProbs)
      #}
      expertWs <- expertNextWs
    }else{
      expertCumulativeLosses <- expertCumulativeLosses + expertLosses
      curExpert <- which.min(expertCumulativeLosses)
    }
    
  }
  return(list(alphaSeq,errSeqAdapt,errSeqFixed,gammaSeq,meanAlphaSeq,meanErrSeq,meanGammas))
}




