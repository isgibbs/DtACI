### Code for reproducing COVID-19 experiments found in https://arxiv.org/abs/2208.08401. The prediction method and associated
### code used here is adapted from work by Ryan Tibshirani https://delphi.cmu.edu/blog/2020/09/21/can-symptoms-surveys-improve-covid-19-forecasts/

library(covidcast)
library(dplyr)
library(tidyr)

################# Construct the Dataset #########################

# Function to append shift values (lags or leads) to data frame
append_shifts = function(df, shifts) {
  # Make sure that we have a complete record of dates for each geo_value (fill
  # with NAs as necessary)
  df_all = df %>% group_by(geo_value) %>%
    summarize(time_value = seq.Date(as.Date(min(time_value)),
                                    as.Date(max(time_value)),
                                    by = "day")) %>% ungroup()
  df = full_join(df, df_all, by = c("geo_value", "time_value"))
  
  # Group by geo value, sort rows by increasing time
  df = df %>% group_by(geo_value) %>% arrange(time_value)
  
  # Load over shifts, and add lag value or lead value
  for (shift in shifts) {
    fun = ifelse(shift < 0, dplyr::lag, dplyr::lead)
    varname = sprintf("value%+d", shift)
    df = mutate(df, !!varname := fun(value, n = abs(shift)))
  }
  # Ungroup and return
  return(ungroup(df))
}

# Some useful functions for transformations
Log = function(x, a = 0.01) log(x + a)
Exp = function(y, a = 0.01) exp(y) - a
Logit = function(x, a = 0.01) log((x + a) / (1 - x + a))
Sigmd = function(y, a = 0.01) (exp(y) * (1 + a) - a) / (1 + exp(y))
Id = function(x) x

#### Parameters #####

# Transforms to consider, in what follows
trans = Id
inv_trans = Id

# Rescale factors for our signals: bring them all down to proportions (between
# 0 and 1)
rescale_g = 1e-2 # Originally a percentage
rescale_f = 1e-2 # Originally a percentage
rescale_c = 1e-5 # Originally a count per 100,000 people

n = 14 # Number of trailing days to use for training set
lp_solver = "glpk" # LP solver to use in quantile_lasso()
verbose = TRUE # Print intermediate progress to console?

#### Data #####
geo_values = covidcast_signal("jhu-csse", "confirmed_cumulative_num",
                              "2020-05-14", "2020-05-14")  %>% pull(geo_value)

# Fetch county-level Facebook % CLI-in-community signals, and JHU
# confirmed case incidence proportion
start_day = "2020-04-11"
end_day = "2022-05-01"
f = covidcast_signal("fb-survey", "smoothed_hh_cmnty_cli",
                     start_day, end_day) %>%
  filter(geo_value %in% geo_values) %>%
  select(geo_value, time_value, value)
c = covidcast_signal("jhu-csse", "confirmed_7dav_incidence_prop",
                     start_day, end_day) %>%
  filter(geo_value %in% geo_values) %>%
  select(geo_value, time_value, value)

# Find "complete" counties, present in all three data signals at all times
geo_values_complete = intersect(f$geo_value,c$geo_value)

# Filter to complete counties, transform the signals, append 1-2 week lags to
# all three, and also 1-2 week leads to case rates
lags = 1 * c(-7,-14)
leads = 1 * 7
fm = f %>% filter(geo_value %in% geo_values_complete) %>%
  mutate(value = trans(value * rescale_f)) %>%
  append_shifts(shifts = lags)
cm = c %>% filter(geo_value %in% geo_values_complete) %>%
  mutate(value = trans(value * rescale_c)) %>%
  append_shifts(shifts = c(lags, leads))

# Rename columns
colnames(fm) = sub("^value", "fb", colnames(fm))
colnames(cm) = sub("^value", "case", colnames(cm))

# Make one big matrix by joining these two data frames
z = full_join(fm,cm, by = c("geo_value", "time_value"))

##### Analysis #####
res_list = vector("list", length = length(leads))

# Loop over lead, forecast dates, build models and record errors (warning: this
# computation takes a while)
for (i in 1:length(leads)) {
  lead = leads[i]; if (verbose) cat("***", lead, "***\n")
  
  # Create a data frame to store our forecast results. Code below populates its
  # rows in a way that breaks from typical dplyr operations, done for efficiency
  res_list[[i]] = z %>%
    filter(between(time_value, as.Date(start_day) - min(lags) + lead,
                   as.Date(end_day) - lead)) %>%
    select(geo_value, time_value) %>%
    mutate(err0 = as.double(NA),err1 = as.double(NA), lead = lead)
  valid_dates = unique(res_list[[i]]$time_value)
  
  for (k in 1:length(valid_dates)) {
    date = valid_dates[k]; if (verbose) cat(format(date), "... ")
    
    # Filter down to training set and test set
    z_tr = z %>% filter(between(time_value, date - lead - n, date - lead))
    z_te = z %>% filter(time_value == date)
    inds = which(res_list[[i]]$time_value == date)
    
    # Create training and test responses
    y_tr = z_tr %>% pull(paste0("case+", lead))
    y_te = z_te %>% pull(paste0("case+", lead))
    
    # Strawman model
    if (verbose) cat("0")
    y_hat = z_te %>% pull(case)
    res_list[[i]][inds,]$err0 = abs(inv_trans(y_hat) - inv_trans(y_te))
    
    # Cases and Facebook model
    x_tr_case = z_tr %>% select(starts_with("case") & !contains("+"))
    x_te_case = z_te %>% select(starts_with("case") & !contains("+"))
    x_tr_fb = z_tr %>% select(starts_with("fb"))
    x_te_fb = z_te %>% select(starts_with("fb"))
    x_tr = cbind(x_tr_case, x_tr_fb)
    x_te = cbind(x_te_case, x_te_fb)
    ok = complete.cases(x_tr, y_tr)
    if (sum(ok) > 0) {
      obj = lm(y_tr[ok] ~ as.matrix(x_tr[ok,]))
      y_hat = as.matrix(x_te)%*%(obj$coefficients[-1]) + obj$coefficients[1]
      res_list[[i]][inds,]$err1 = as.numeric(abs(inv_trans(y_hat) - inv_trans(y_te)))
    }
  }
}

# Bind results over different leads into one big data frame, and save
res = do.call(rbind, res_list)
#save(list = ls(), file = "demo.rda")

########################## Run FACI #################################

#load("forecast-demo/demo.rda")
library(dplyr)
library(tidyr)
library(ggplot2)
source("FACI.R")

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

### Compute beta_ts for the given conformity scores. Lookback specifies the amount of time periods used to construct the model,
### epsilon specifies a numerical error tolerence for the binary search for beta_t,
### geosToUse specifies which counties we compute the betas for (NULL means look at all counties)
computeBetasByGeoByTime <- function(scores,lookback=1,epsilon=0.001,geosToUse = NULL){
  dates = unique(scores$time_value)
  geovals = unique(scores$geo_value)
  T <- length(dates)
  if(is.null(geosToUse)){
    geosToUse <- geovals
  }
  
  ### Initialize data storage variable
  betaSeqMat <- matrix(0,nrow = length(geosToUse),ncol=length(dates)-1)
  
  for(t in 2:length(dates)){
    prevScores <- scores$err1[scores$time_value==dates[t-1]]
    for(i in 1:length(geosToUse)){
      newScore <- scores$err1[scores$time_value==dates[t] & scores$geo_value==geosToUse[i]]
      if(length(newScore) > 0){
        betaSeqMat[i,t-1] <- findBeta(prevScores,newScore,epsilon)
      }else{
        betaSeqMat[i,t-1] <- NA
      }
    }
    if(t %% 10 == 0){
      print(sprintf("Done %g steps",t))
    }
  }
  
  return(betaSeqMat)
}

# Calculate the scaled errors relative to the strawman's error
res_final <- res %>%
  drop_na() %>%                                       # Restrict to common time
  mutate(err1 = err1 / err0) %>% 
  ungroup() %>%
  select(-err0)
res_final <- res_final %>% group_by(time_value) %>% arrange(geo_value) 

### Compute all beta_t values
alpha <- 0.1
gammaGrid <- c(0.001,0.002,0.004,0.008,0.0160,0.032,0.064,0.128)
geosToUse <- c("06075","36061","12086", "48113")  ## SF, NY, Miami-Dade, Dallas

allBetas <- computeBetasByGeoByTime(res_final,geosToUse = geosToUse)
rownames(allBetas) <- c("06075","36061","12086", "48113")

### Compute and plot final results
library(latex2exp)
library(cowplot)
library(grid)
source("PlottingCode.R")

myDates <- unique(res_final$time_value)[-1]
countiesToUse <- c("06075","36061","12086", "48113")  ## SF, NY, Miami-Dade, Dallas

allRes <- list()
myPlots <- list()
minRangeY <- 1
maxRangeY <- 0
count <- 1
for(geo in countiesToUse){
  allRes[[count]] <- conformalAdaptStable(allBetas[rownames(allBetas) ==geo],alpha,gammaGrid,sigma=1/500,eta=2.8)
  myPlots[[count]] <- plotLocalAvgDatesV2(allRes[[count]][[5]],allRes[[count]][[6]],allRes[[count]][[3]],alpha,myDates,startUp=100) 
  if(ggplot_build(myPlots[[count]][[2]] )$layout$panel_params[[1]]$y.range[2] > maxRangeY){
    maxRangeY <- ggplot_build(myPlots[[count]][[2]] )$layout$panel_params[[1]]$y.range[2]
  }
  if(ggplot_build(myPlots[[count]][[2]] )$layout$panel_params[[1]]$y.range[1] < minRangeY){
    minRangeY <- ggplot_build(myPlots[[count]][[2]] )$layout$panel_params[[1]]$y.range[1]
  }
  count <- count + 1
}
maxRangeY<-maxRangeY+0.05
for(i in 1:length(countiesToUse)){
  myPlots[[i]][[2]] <- myPlots[[i]][[2]] + ylim(minRangeY,maxRangeY) + theme(axis.text=element_text(size=11.5),axis.title=element_text(size=15,face="plain"))
}

legend <- get_legend(myPlots[[1]][[2]]+theme(legend.position = "top",legend.text=element_text(size=15),legend.title=element_blank()))
blankplot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()
myPlots[[1]][[2]] <- myPlots[[1]][[2]] +  annotate('text',x=as.Date("2020-10-22"), y=maxRangeY, label="San Francisco, CA",size=6,fontface="bold") + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) + theme(legend.position = "none")
myPlots[[2]][[2]] <- myPlots[[2]][[2]]  +  annotate('text',x=as.Date("2020-10-05"), y=maxRangeY, label="New York, NY",size=6,fontface="bold") + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())  + theme(legend.position = "none")
myPlots[[3]][[2]] <- myPlots[[3]][[2]] +  annotate('text',x=as.Date("2020-10-08"), y=maxRangeY, label="Miami-Dade, FL",size=6,fontface="bold") +  theme(axis.title.x = element_blank()) + theme(legend.position = "none")
myPlots[[4]][[2]] <- myPlots[[4]][[2]] +  annotate('text',x=as.Date("2020-09-16"), y=maxRangeY, label="Dallas, TX",size=6,fontface="bold") + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()) +  theme(axis.title.x = element_blank())  + theme(legend.position = "none")

grid.arrange(legend,blankplot,myPlots[[1]][[2]],myPlots[[2]][[2]],myPlots[[3]][[2]],myPlots[[4]][[2]],ncol=2,nrow=3,widths=c(3.5,3) , heights = c(0.2, 2.5,2.5),bottom=textGrob("Time", gp=gpar(fontsize=20)))



