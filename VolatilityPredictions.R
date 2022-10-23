### This file contains code to reproduce the volatility prediction experiments found in https://arxiv.org/abs/2208.08401 

source("FACI.R")
library(rugarch)

### Compute sequence of conformity scores by fitting Garch(garchP,garchQ) model. 
### Returns is a sequence of stock market returns, lookback specifies the dataset \{R_s\}_{t-lookback+1 \leq s \leq t-1}
### used to fit the model and make predictions at time t, startup specifies the first time we make a prediction set,
### badscores = FALSE means we use the conformity scores S_t := |V_t- \hat{\sigma}^2_t|/\hat{\sigma}^2_t else we use
### S_t := |V_t- \hat{\sigma}^2_t|
garchConformalForcastingComputeScores <- function(returns,lookback=1250,garchP=1,garchQ=1,startUp = 100,badScores=FALSE){
  T <- length(returns)
  startUp <- max(startUp,lookback)
  garchSpec <- ugarchspec(mean.model=list(armaOrder = c(0, 0),include.mean=FALSE),variance.model=list(model="sGARCH",garchOrder=c(1,1)),distribution.model="norm")
  scores <- rep(0,T-startUp+1)
  sigmaSeq <- rep(0,T-startUp+1)
  
  for(t in startUp:T){
    garchFit <- ugarchfit(garchSpec, returns[(t-lookback+1):(t-1) ],solver="hybrid")
    sigmaNext <- sigma(ugarchforecast(garchFit,n.ahead=1))
    if(!badScores){
      scores[t-startUp + 1] <- abs(returns[t]^2- sigmaNext^2)/sigmaNext^2
    }else{
      scores[t-startUp + 1] <- abs(returns[t]^2- sigmaNext^2)
    }
    sigmaSeq[t-startUp + 1] <- sigmaNext^2
    
    if(t %% 100 == 0){
      print(sprintf("Done %g steps",t))
    }
  }
  return(list(scores,sigmaSeq))
}

### Use a binary search to find the lowest quantile of recentScores that is above curScore
### Epsilon gives numerical error tolerance in the binary search
findBeta <- function(recentScores,curScore,epsilon=0.001){
  top <- 1
  bot <- 0
  mid <- (top+bot)/2
  while(top-bot > epsilon){
    if(quantile(recentScores,1-mid)>curScore){
      bot <- mid
      mid <- (top+bot)/2
    }else{
      top <- mid
      mid <- (top+bot)/2
    }
  }
  return(mid)
}

### Given a sequence of conformity scores compute the corresponding values for beta_t
### Epsilon gives numerical error tolerance in a binary search
garchConformalForcastingComputeBetas <- function(scores,lookback=Inf,epsilon=0.001){
  T <- length(scores)
  ### Initialize data storage variables
  betaSeq <- rep(0,T-1)
  
  for(t in 2:T){
    recentScores <- scores[max(t- lookback,1):(t-1)]
    betaSeq[t-1] <- findBeta(recentScores,scores[t],epsilon)
    if(t %% 1000 == 0){
      print(sprintf("Done %g steps",t))
    }
  }
  
  return(betaSeq)
}

### Get stock market price data
amd <- read.csv("MajorDailyReturns/AMDDailyPrices.csv",header=TRUE)
nividia <- read.csv("MajorDailyReturns/NividiaDailyPrices.csv",header=TRUE)
blackberry <- read.csv("MajorDailyReturns/BlackBerryDailyPrices.csv",header=TRUE)
fannie <- read.csv("MajorDailyReturns/FannieMaeDailyPrices.csv",header=TRUE)
pricesAMD <- rev(amd$Open)
returnsAMD <- sapply(2:length(pricesAMD),function(x){(pricesAMD[x] - pricesAMD[x-1])/pricesAMD[x-1]})
returnsAMD <- returnsAMD[(length(returnsAMD)-7500):length(returnsAMD)]
pricesBlackberry<- rev(blackberry$Open)
returnsBlackberry <- sapply(2:length(pricesBlackberry),function(x){(pricesBlackberry[x] - pricesBlackberry[x-1])/pricesBlackberry[x-1]})
pricesNividia <- rev(nividia$Open)
returnsNividia <- sapply(2:length(pricesNividia),function(x){(pricesNividia[x] - pricesNividia[x-1])/pricesNividia[x-1]})
pricesFannie <- rev(fannie$Open)
returnsFannie <- sapply(2:length(pricesFannie),function(x){(pricesFannie[x] - pricesFannie[x-1])/pricesFannie[x-1]})
returnsFannie <- returnsFannie[(length(returnsFannie)-7500):length(returnsFannie)]

### Compute all conformity scores
lookback <- 1250
amdScores <- garchConformalForcastingComputeScores(returnsAMD,lookback=lookback)
blackberryScores <- garchConformalForcastingComputeScores(returnsBlackberry,lookback=lookback)
nividiaScores <- garchConformalForcastingComputeScores(returnsNividia,lookback=lookback)
fannieScores <- garchConformalForcastingComputeScores(returnsFannie,lookback=lookback)

amdBadScores <- garchConformalForcastingComputeScores(returnsAMD,lookback=lookback,badScores=TRUE)
blackberryBadScores <- garchConformalForcastingComputeScores(returnsBlackberry,lookback=lookback,badScores=TRUE)
nividiaBadScores <- garchConformalForcastingComputeScores(returnsNividia,lookback=lookback,badScores=TRUE)
fannieBadScores <- garchConformalForcastingComputeScores(returnsFannie,lookback=lookback,badScores=TRUE)

### Compute betas
amdBetas <- garchConformalForcastingComputeBetas(amdScores[[1]],lookback=lookback)
blackberryBetas <- garchConformalForcastingComputeBetas(blackberryScores[[1]],lookback=lookback)
nividiaBetas <- garchConformalForcastingComputeBetas(nividiaScores[[1]],lookback=lookback)
fannieBetas <- garchConformalForcastingComputeBetas(fannieScores[[1]],lookback=lookback)

amdBadBetas <- garchConformalForcastingComputeBetas(amdBadScores[[1]],lookback=lookback)
blackberryBadBetas <- garchConformalForcastingComputeBetas(blackberryBadScores[[1]],lookback=lookback)
nividiaBadBetas <- garchConformalForcastingComputeBetas(nividiaBadScores[[1]],lookback=lookback)
fannieBadBetas <- garchConformalForcastingComputeBetas(fannieBadScores[[1]],lookback=lookback)


### Run FACI
gammaGrid <- c(0.001,0.002,0.004,0.008,0.0160,0.032,0.064,0.128)
amdStableGridGammas <- conformalAdaptStable(amdBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
blackberryStableGridGammas <- conformalAdaptStable(blackberryBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
nividiaStableGridGammas <- conformalAdaptStable(nividiaBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
fannieStableGridGammas <- conformalAdaptStable(fannieBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)

amdStableGridGammasBad <- conformalAdaptStable(amdBadBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
blackberryStableGridGammasBad <- conformalAdaptStable(blackberryBadBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
nividiaStableGridGammasBad <- conformalAdaptStable(nividiaBadBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)
fannieStableGridGammasBad <- conformalAdaptStable(fannieBadBetas,alpha,gammaGrid,sigma=1/500,eta=2.8)

### Plot some of the results
source("PlottingCode.R")

allDates <- list(amd$Date[1:length(amdStableGridGammasBad[[1]])])
for(d in 1:length(allDates)){
  allDates[[d]] <- rev(allDates[[d]])
  allDates[[d]] <- as.Date(allDates[[d]], '%m/%d/%Y')
}

for(d in 1:length(allDates)){
  allDates[[d]] <- gsub("009","199",allDates[[d]])
  allDates[[d]] <- gsub("000","200",allDates[[d]])
  allDates[[d]] <- gsub("001","201",allDates[[d]])
  allDates[[d]] <- gsub("002","202",allDates[[d]])
  allDates[[d]] <- gsub("2201","2001",allDates[[d]])
  allDates[[d]] <- gsub("2202","2002",allDates[[d]])
  allDates[[d]] <- gsub("0199","2009",allDates[[d]])
}
myplots <- plotLocalAvgDatesV2(amdResults[[1]],amdResults[[2]],amdResults[[3]],alpha,allDates[[1]],bandwidth=1,startUp=250)
grid.arrange(myplots[[1]],myplots[[2]],ncol=2)

