# This script runs the subperiod forward selection model at eight-weeek horizon.
# Also it saves all the regression results in a single word file for Table A.VII. 
# Output
# forwardSelection_subperiod_all_8wks.doc

# ------------------------------------------------------------------------
# Clear Environment 
rm(list = ls())

# Libraries 
library(sjPlot)
library(lmtest)
library(dplyr)
library(tidyr)
library(sandwich)
library(readstata13)
library(selectiveInference)
library(broom)
library(purrr)
library(data.table)
library(gtools)

# source("/if/appl/R/Functions/IFfunctions.r") # loads necessary packages and IF helper functions

# Set working directory
dir <- "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/"
setwd(dir)

# Subperiod 1 - 9


#### sub1

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("1998-04-01") & date_Fri <= as.Date("1999-11-30"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("1998-04-01") & date_Tue <= as.Date("1999-11-30"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Forward selection function                                                                                     
# ------------------------------------------------------------------------

fwdSelection <- function(subset, depVar, horizon){
  
  results <- list()
  
  # Regress all variables on lagged dependent variables 
  detrended.lagDepVar.adjusted <- subset
  for (i in 1:ncol(subset)) {
    if(i == 1){
      lm <- lm(subset[,1] ~ dplyr::lag(subset[,2], horizon), data = subset, na.action=na.exclude)
      pred.value <- predict.lm(lm)
      detrended.lagDepVar.adjusted[,1] <- subset[,1] - pred.value
    } else {
      lm <- lm(subset[,i] ~ subset[,2], data = subset, na.action=na.exclude)
      pred.value <- predict.lm(lm)
      detrended.lagDepVar.adjusted[,i] <- subset[,i] - pred.value
    }
  }
  
  # Remove lagged dependent variable
  detrended.lagDepVar.adjusted <- dplyr::select(detrended.lagDepVar.adjusted, -c(2))
  
  # Lag detrended and lagged-dependent-variable adjusted independent variables for forward selection
  final.df <- detrended.lagDepVar.adjusted
  for(i in 2:ncol(detrended.lagDepVar.adjusted)){
    var <- colnames(detrended.lagDepVar.adjusted)[i]
    if(sum(var == toLag) == 1){
      final.df[,i] <- dplyr::lag(detrended.lagDepVar.adjusted[,i], horizon)
    }
  }
  
  # Remove risk premia measures from final.df
#  final.df.less.risk <- dplyr::select(final.df, -c("ovx_diff", "sdf_fullSample"))
  
  subset <- as.matrix(final.df)
  
  # Omit NAs from detrended matrix
  subset <- na.omit(subset)
  
  # Select x matrix and y vector of variables 
  x <- as.matrix(subset[,2:ncol(subset)])
  y <- as.vector(subset[,1])
  
  # Run forward stepwise, plot results
  fsfit <- fs(x,y,7)
  
  # compute sequential p-values and confidence intervals
  output <- fsInf(fsfit)
  
  index <- output$vars
  selected <- colnames(x)[index]
  
  lm.df <- as.data.frame(cbind(y,x[,index]))
  
  results$lm.0 <- lm(y ~ ., data = lm.df)
  results$sd <- c(sd(lm.df[,1]))
  results$sd <- c(results$sd, c(sd(lm.df[,2])))
  results$sd <- c(results$sd, c(sd(lm.df[,3])))
  results$sd <- c(results$sd, c(sd(lm.df[,4])))
  results$sd <- c(results$sd, c(sd(lm.df[,5])))
  results$sd <- c(results$sd, c(sd(lm.df[,6])))
  results$sd <- c(results$sd, c(sd(lm.df[,7])))
  results$sd <- c(results$sd, c(sd(lm.df[,8])))  
  results$stdcoef <- c(results$lm.0$coefficients[1],results$lm.0$coefficients[2:8]/results$sd[1]*results$sd[-1])
  # There is discrepancy between p-values of standardized regressions and plain regressions due to the centering and scaling of data to have unit variance
  # Since we want to report standardized coefficients with p-values from original regressions, 
  # we create a customized table which cannot be created by tab_model or other R functions.
  # tab_model automatically calculates p-values of the regression coefficients,
  # so we cannot customize within tab_model to report std. coefficients with the *, **, *** from original regression coefficients.
  # tab_model argument "show.std= T" shows the std. coefficients in a separate column. But we cannot put original regression stars(*) next to it. 
  
  #Create a function to convert p-value to stars
  makeStars <- function(x){
    stars <- c("***", "**", "*", "")
    vec <- c(0, 0.01, 0.05, 0.1)
    i <- findInterval(x, vec)
    stars[i]
  }
  
  
  results$coeff_pval = data.table(Predictors=names(results$lm.0$coefficients), Estimates=paste0(round(results$stdcoef, digits=2)," ",
                                                                                               makeStars(summary(results$lm.0)$coefficients[,4])))
  results$coeff_pval <- rbind(results$coeff_pval,list("Observations", length(results$lm.0$model$y)))
  results$coeff_pval <- rbind(results$coeff_pval,list("Rsq/Adj.Rsq", paste0(round(summary(results$lm.0)$r.squared,2), "/",
                                                            round(summary(results$lm.0)$adj.r.squared,2))))
  results$coeff_pval <-results$coeff_pval[-1,] #erase the row for intercept
  
  

  return(results)
}

# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub1 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub1 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub1 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub1 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub1 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub1 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub1 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub1 <- fwdSelection(subset, depVar, horizon)



#### sub2

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("1999-12-01") & date_Fri <= as.Date("2002-07-31"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("1999-12-01") & date_Tue <= as.Date("2002-07-31"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub2 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub2 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub2 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub2 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub2 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub2 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub2 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub2 <- fwdSelection(subset, depVar, horizon)


#### sub3

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2002-08-01") & date_Fri <= as.Date("2005-03-31"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2002-08-01") & date_Tue <= as.Date("2005-03-31"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub3 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub3 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub3 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub3 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub3 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub3 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub3 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub3 <- fwdSelection(subset, depVar, horizon)


#### sub4

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2005-04-01") & date_Fri <= as.Date("2007-11-30"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2005-04-01") & date_Tue <= as.Date("2007-11-30"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub4 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub4 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub4 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub4 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub4 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub4 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub4 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub4 <- fwdSelection(subset, depVar, horizon)


#### sub5

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2007-12-01") & date_Fri <= as.Date("2009-06-30"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2007-12-01") & date_Tue <= as.Date("2009-06-30"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub5 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub5 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub5 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub5 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub5 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub5 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub5 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub5 <- fwdSelection(subset, depVar, horizon)


#### sub6

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2009-07-01") & date_Fri <= as.Date("2012-02-29"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2009-07-01") & date_Tue <= as.Date("2012-02-29"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub6 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub6 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub6 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub6 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub6 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub6 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub6 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub6 <- fwdSelection(subset, depVar, horizon)



#### sub7

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2012-03-01") & date_Fri <= as.Date("2014-10-31"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2012-03-01") & date_Tue <= as.Date("2014-10-31"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub7 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub7 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub7 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub7 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub7 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub7 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub7 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub7 <- fwdSelection(subset, depVar, horizon)


#### sub8

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2014-11-01") & date_Fri <= as.Date("2017-06-30"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2014-11-01") & date_Tue <= as.Date("2017-06-30"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub8 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub8 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub8 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub8 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub8 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub8 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub8 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub8 <- fwdSelection(subset, depVar, horizon)


#### sub9

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("2017-07-01") & date_Fri <= as.Date("2020-03-31"))


# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.prices <- subset.prices
for (i in 1:ncol(subset.prices)) {
  lm <- lm(subset.prices[,i] ~ trend, data = subset.prices, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.prices[,i] <- subset.prices[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.prices <- dplyr::select(detrended.prices, -c(trend))


# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX
data.physical <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta") 
data.physical <- subset(data.physical, date_Tue>= as.Date("2017-07-01") & date_Tue <= as.Date("2020-03-31"))

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))


# ------------------------------------------------------------------------
# Define list of variables we want to lag in MC.FS function                                                                                   
# ------------------------------------------------------------------------
toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX",
           "vix_diff", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")




# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

FutRet_t8_sub9 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

DSpot_t8_sub9 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil Volatility Changes - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, 
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

DOilVol_t8_sub9 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

xomRet_t8_sub9 <- fwdSelection(subset, depVar, horizon)


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

bpRet_t8_sub9 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

rdsaRet_t8_sub9 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

DInv_t8_sub9 <- fwdSelection(subset, depVar, horizon)

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

DProd_t8_sub9 <- fwdSelection(subset, depVar, horizon)







#############################################################################
# ------------------------------------------------------------------------
#  Output subperiod regression results using stepwise forward selection                                                                              
# ------------------------------------------------------------------------
# We record all results in a single document for convenience. 
library(rtf)
rtffile <- RTF("20210821/forwardSelection_subperiod_all_8wks.doc") 

#FutRet
addHeader(rtffile,title = "Table A.VII", subtitle ="Panel VII.1 : Predicting 8-weeks ahead FutRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#DSpot
addHeader(rtffile,title = "", subtitle ="Panel VII.2 : Predicting 8-weeks ahead DSpot in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#DOilVol
addHeader(rtffile,title = "", subtitle ="Panel VII.3 : Predicting 8-weeks ahead DOilVol in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#xomRet
addHeader(rtffile,title = "", subtitle ="Panel VII.4 : Predicting 8-weeks ahead xomRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#bpRet
addHeader(rtffile,title = "", subtitle ="Panel VII.5 : Predicting 8-weeks ahead bpRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#rdsaRet
addHeader(rtffile,title = "", subtitle ="Panel VII.6 : Predicting 8-weeks ahead rdsaRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#DInv
addHeader(rtffile,title = "", subtitle ="Panel VII.7 : Predicting 8-weeks ahead DInv in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

#DProd
addHeader(rtffile,title = "", subtitle ="Panel VII.8 : Predicting 8-weeks ahead DProd in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_t8_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "* p<0.1 ** p<0.05 *** p<0.01", bold=FALSE, italic=TRUE)
addText(rtffile, "\n")

done(rtffile)


















