### Code for reproducing block bootstrap error bars found in https://arxiv.org/abs/2208.08401

source("VolatilityPredictions.R")

###################### Definitions of Key Functions ####################

source("VolatilityPredictions.R")

### Run one block bootstrap iteration. Inputs are blocks of the data, sequence of stock returns, amount of time to lookback
### when constructing the model, target level alpha, size of each block, grid of gamma values for FACI,
### bands B_i in which we should compute the average coverage, width of the bands, and a flag indicating whether
### we should use the conformity score S_t := |V_t- \hat{\sigma}^2_t|/\hat{\sigma}^2_t (goodScores = TRUE) or 
### S_t := |V_t- \hat{\sigma}^2_t| (goodScores = FALSE)
### Return value is the observed average coverage in each band
oneBootRun <- function(blocks,returns,lookback,alpha,blockSize,gammaGrid,bands,bandwidth,goodScores=TRUE){
    blocksBoot <- sample(blocks,size = length(blocks),replace=TRUE) 
    indexSeq <- as.vector(sapply(blocksBoot,function(x){x:(x+blockSize-1)}))
    scores <- garchConformalForcastingComputeScores(returns[indexSeq],lookback=lookback,badScores=!goodScores)
    betas <- garchConformalForcastingComputeBetas(scores[[1]],lookback=lookback)
    res <- conformalAdaptStable(betas,alpha,gammaGrid,sigma=1/1000,eta=2.72)
    return(1-sapply(bands,function(x){mean(res[[6]][res[[5]]> x & res[[5]] < x+bandwidth])}))
}

### Run the block bootstrap to compute error bars. First four inputs are the observed sequences alphabar_t, 
### err_t(alphabar_t) (i.e. the errors obtained with alphabar_t), the sequence of stock market returns, 
### and the grid of gamma values used in FACI. nBootSamples is the total number of bootstrap samples to compute. 
### Descriptions of all  other inputs can be found in the previous comment.
bootForErrs <- function(alphats,errts,returns,gammaGrid = c(0.001,0.002,0.004,0.008,0.0160,0.032,0.064,0.128),
                        lookback = 1250,bandwidth=0.02,alpha=0.1,bands=NULL,blockSize = 100,
                        nBootSamples=100,goodScores=TRUE){
  if(is.null(bands)){
    bands = seq(min(alphats),max(alphats),by=bandwidth)
  }
  
  returns <- returns[1:(length(returns) -  ( length(returns) %% blockSize))]
  blocks <- seq(1,length(returns),by=blockSize)
  
  allFreqs <- matrix(0,nrow = length(bands),ncol = nBootSamples)
  for(i in 1:nBootSamples){
    allFreqs[,i] <- oneBootRun(blocks,returns,lookback,alpha,blockSize,gammaGrid,bands,bandwidth,goodScores)
  }
  return(allFreqs)
}

### Compute quantiles of the bootstrap distribution allFreqs output by bootForErrs.
computeBootstrapQuantiles <- function(alphats,errts,bands,bandwidth,allFreqs,coverageProb=0.1){
  coverageFrequencies <- 1-sapply(bands,function(x){mean(errts[alphats> x & alphats < x+bandwidth])})
  return(sapply(1:nrow(allFreqs),function(i){c(quantile(allFreqs[i,], coverageProb/2,na.rm=TRUE), 
                                               quantile(allFreqs[i,], 1-coverageProb/2,na.rm=TRUE))}))
}

### Code for plotting the final results. Inputs are alphabar_t, err_t(alphabar_t) (i.e. the errors obtained with alphabar_t), 
### the quantiles output by computeBootstrapQuantiles and cutoff value ymin that limits the range of the y-axis to [ymin,1].
plotCoverageByAlphat <- function(alphats,errts,errorQuants,ymin=0.5){
  coverageSampleSizes <- sapply(bands,function(x){length(alphats[alphats> x & alphats < x+bandwidth])})
  coverageFrequencies <- 1-sapply(bands,function(x){mean(errts[alphats> x & alphats < x+bandwidth])})
  midpoints <- bands + bandwidth/2
  
  myData <- data.frame(alphat=midpoints,coverageFreq = coverageFrequencies,bot = errorQuants[1,], top =  errorQuants[2,])
  myData$top[is.na(coverageFrequencies)] = NA
  myData$bot[is.na(coverageFrequencies)] = NA
  myData$top <- sapply(myData$top, function(x){min(x,1)})
  myData$bot <- sapply(myData$bot, function(x){max(x,ymin)})
  myData[coverageSampleSizes<25,] <- c(NA,NA,NA,NA)
  print(myData)
  
  p <- ggplot(myData,aes(alphat,coverageFreq)) + geom_point(aes(alphat,coverageFreq),color="red") +   geom_errorbar(aes(ymin=bot, ymax=top), width=bandwidth/2) 
  p <- p + geom_hline(yintercept=1-alpha, linetype="dashed", color = "black",size = 1.2) +
    ylab("Empirical Conditional Coverage") + xlab(TeX("$\\bar{\\alpha}_t$"))+
    theme(axis.text=element_text(size=11.5),axis.title=element_text(size=15,face="plain"))
  
  print(coverageSampleSizes)
  return(p)
}

###################### Run the bootstrap and plot results ####################
### The data here comes from the output of the file VolatilityPredictions.R

### Compute partition of [0,1] within which we will compute the conditional coverages
minAlphat <- min(c(amdStableGridGammas[[5]],blackberryStableGridGammas[[5]],nividiaStableGridGammas[[5]],fannieStableGridGammas[[5]]))
maxAlphat <- max(c(amdStableGridGammas[[5]],blackberryStableGridGammas[[5]],nividiaStableGridGammas[[5]],fannieStableGridGammas[[5]]))
bandwidth = 0.02
bands = seq(minAlphat,maxAlphat,by=bandwidth)

### Compute all error bars and plot results
amdAllFreqs <- bootForErrs(amdStableGridGammas[[5]],amdStableGridGammas[[6]],returnsAMD,bands=bands,bandwidth=bandwidth)
amdErrorQuants <- computeBootstrapQuantiles(amdStableGridGammas[[5]],amdStableGridGammas[[6]],bands,bandwidth,amdAllFreqs)
amdAlphaVersusCoveragePlot <- plotCoverageByAlphat(amdStableGridGammas[[5]],amdStableGridGammas[[6]],amdErrorQuants) +
  annotate('text',x=minAlphat+0.0132, y=1.25, label="AMD",size=6,fontface="bold") + ylim(0.5,1.25) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())

bbAllFreqs <- bootForErrs(blackberryStableGridGammas[[5]],blackberryStableGridGammas[[6]],returnsBlackberry,
                          bands=bands,bandwidth=bandwidth)
bbErrorQuants <- computeBootstrapQuantiles(blackberryStableGridGammas[[5]],blackberryStableGridGammas[[6]],
                                           bands,bandwidth,bbAllFreqs)
bbAlphaVersusCoveragePlot <- plotCoverageByAlphat(blackberryStableGridGammas[[5]],
                                                  blackberryStableGridGammas[[6]],bbErrorQuants) +
  annotate('text',x=minAlphat+0.015, y=1.25, label="BlackBerry",size=6,fontface="bold") + ylim(0.5,1.25)

nividiaAllFreqs <- bootForErrs(nividiaStableGridGammas[[5]],nividiaStableGridGammas[[6]],returnsNividia,
                               bands=bands,bandwidth=bandwidth)
nividiaErrorQuants <- computeBootstrapQuantiles(nividiaStableGridGammas[[5]],nividiaStableGridGammas[[6]],
                                                bands,bandwidth,nividiaAllFreqs)
nividiaAlphaVersusCoveragePlot <- plotCoverageByAlphat(nividiaStableGridGammas[[5]],nividiaStableGridGammas[[6]],
                                                       nividiaErrorQuants) +
  annotate('text',x=minAlphat+0.028, y=1.25, label="Nvidia",size=6,fontface="bold") + ylim(0.5,1.25) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())

fannieAllFreqs <- bootForErrs(fannieStableGridGammas[[5]],fannieStableGridGammas[[6]],returnsFannie,
                              bands=bands,bandwidth=bandwidth)
fannieErrorQuants <- computeBootstrapQuantiles(fannieStableGridGammas[[5]],fannieStableGridGammas[[6]],
                                               bands,bandwidth,fannieAllFreqs)
fannieAlphaVersusCoveragePlot <- plotCoverageByAlphat(fannieStableGridGammas[[5]],fannieStableGridGammas[[6]],fannieErrorQuants) +
  annotate('text',x=minAlphat+0.016, y=1.25, label="Fannie Mae",size=6,fontface="bold") + ylim(0.5,1.25)+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())


grid.arrange(nividiaAlphaVersusCoveragePlot,amdAlphaVersusCoveragePlot,bbAlphaVersusCoveragePlot,fannieAlphaVersusCoveragePlot,ncol=2)


