library(ggplot2)
library(gridExtra)

set.seed(12345)

### Compute an exponentially weighted mean (with weight bandwidth) of elements (index-bottom+1):index of seq.
localMean <- function(seq,index,bandwidth,bottom=300){
  weights <- (bandwidth)^((index-bottom+1):index)
  weights <- weights/sum(weights)
  return(sum(weights*seq[(index-bottom+1):index]))
}

### Inputs are: alphaSequence: the sequence of alphat values, errSeqOC: the sequence of errors produced by FACI, 
### errSeqNC: the sequence of errors produces by alphat = alpha, alpha: the target coverage level, 
### mydates: a vector of dates used to construct the x-axis, 
### bandwidth: A parameter used to compute an exponentially weighted mean, we usually take bandwidth = 1,
### startup: the amount we should look ahead and back in time when computing the local mean.
plotLocalAvgDatesV2 <- function(alphaSequence,errSeqOC,errSeqNC,alpha,mydates,bandwidth=0.995,startUp=300){
  n <- length(errSeqOC)-startUp
  if(!is.null(mydates)){
    temp <- data.frame(alpha = alphaSequence[startUp:n],dates=as.Date(mydates[startUp:n]))
  }else{
    temp <- data.frame(alpha = alphaSequence[startUp:n],dates=startUp:n)
  }
  p1 <- ggplot(temp, aes(dates)) + 
    geom_line(aes(y=alpha, colour="black")) + 
    scale_colour_manual(name = "Method", values =c("black"="black"), labels = c("")) +
    xlab("Time") + ylab("alpha_t") 
  
  covCurveAdapt <- 1-sapply(startUp:n,function(x){localMean(errSeqOC,x+startUp,bandwidth,bottom=2*startUp)})
  covCurveNaive <- 1-sapply(startUp:n,function(x){localMean(errSeqNC,x+startUp,bandwidth,bottom=2*startUp)})
  
  binomseq <- rbinom(length(errSeqOC),1,alpha)
  covCurveBinom <- 1-sapply(startUp:n,function(x){localMean(binomseq,x+startUp,bandwidth,bottom=2*startUp)})
  
  if(!is.null(mydates)){
    temp2 <- data.frame(OnlineUpdate = covCurveAdapt, NaiveMethod = covCurveNaive,binomCurve=covCurveBinom,dates=as.Date(mydates[startUp:n]))
  }else{
    temp2 <- data.frame(OnlineUpdate = covCurveAdapt, NaiveMethod = covCurveNaive,binomCurve=covCurveBinom,dates=startUp:n)
  }
  p2 <- ggplot(temp2,  aes(dates)) + 
    geom_line(aes(y=OnlineUpdate, colour="blue")) + 
    geom_line(aes(y=NaiveMethod, colour="red")) + 
    geom_line(aes(y=binomCurve, colour="z")) +
    scale_colour_manual(name = "Method", values =c("blue"="blue","red"="red","z"="grey"), labels = c("Adaptive Alpha","Fixed Alpha","Bernoulli Sequence")) +
    xlab("Time") + ylab("Local Coverage Level")  + 
    geom_hline(yintercept=1-alpha, linetype="dashed", color = "black",size = 1.2) +
    geom_hline(yintercept=mean(1-errSeqOC), linetype="dashed", color = "blue",size = 1.2) +
    geom_hline(yintercept=mean(1-errSeqNC), linetype="dashed", color = "red",size = 1.2)
  
  return(list(p1,p2))
}
