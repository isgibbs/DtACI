### Determine value of Keps parameter
computeKeps <- function(epsilon=0.01, numerr = 0.0000001){
  total = 0
  i = 0
  while(TRUE){
    newTerm =  1/((i+1)*log(i+2)^(1+epsilon))
    total = total + newTerm
    if(newTerm < numerr){
      return(total)
    }
    i = i +1
  }
}

MVP <- function(scores,alpha = 0.1,m = 40,r = 800000,epsilon = 1, Keps = -1){
  if(Keps == -1){
    Keps = computeKeps(epsilon)
  }
  eta = sqrt(log(m)/(2*Keps*m))
  
  f <- function(n){
    sqrt((n+1)*log(n+2,2)^(1+epsilon))
  }
  
  ns = rep(0,m)
  Vs = rep(0,m)
  Cs = rep(0,m)
  errSeq = rep(0,length(scores))
  alphaSeq = rep(0,length(scores))
  for(i in 1:length(scores)){
    foundFlag = FALSE
    for(j in 1:(m-1)){
      if(Vs[j]*Vs[j+1] <= 0 & !foundFlag){
        foundFlag = TRUE
        if(abs(Cs[j+1]) + abs(Cs[j]) == 0){
          p = 1
        }else{
          p = abs(Cs[j+1])/(abs(Cs[j+1]) + abs(Cs[j]))
        }
        Z = rbinom(1,1,p)
        if(Z == 1){
          threshold = (j)/m - 1/(r*m)
          errSeq[i] = scores[i] > threshold
          ns[j] = ns[j] + 1
          Vs[j] = Vs[j] + (1-errSeq[i] - (1-alpha))
          Cs[j] = (exp(eta*Vs[j]/f(ns[j])) - exp(-eta*Vs[j]/f(ns[j])))/f(ns[j])
        }else{
          threshold = (j)/m 
          errSeq[i] = scores[i] > threshold
          ns[j+1] = ns[j+1] + 1
          Vs[j+1] = Vs[j+1] + (1-errSeq[i] - (1-alpha))
          Cs[j+1] = (exp(eta*Vs[j+1]/f(ns[j+1])) - exp(-eta*Vs[j+1]/f(ns[j+1])))/f(ns[j+1])
        }
      }
    }
    if(!foundFlag){
      if(Vs[1] < 0){
        threshold = 1
        errSeq[i] = 0
        ns[m] = ns[m] + 1
        Vs[m] = Vs[m] + (1-errSeq[i] - (1-alpha))
        Cs[m] = (exp(eta*Vs[m]/f(ns[m])) - exp(-eta*Vs[m]/f(ns[m])))/f(ns[m])
      }else{
        threshold = 0
        errSeq[i] = 1
        ns[1] = ns[1] + 1
        Vs[1] = Vs[1] + (1-errSeq[i] - (1-alpha))
        Cs[1] = (exp(eta*Vs[1]/f(ns[1])) - exp(-eta*Vs[1]/f(ns[1])))/f(ns[1])
      }
    }
    alphaSeq[i] = threshold
  }
  return(list(alphaSeq,errSeq))
}
