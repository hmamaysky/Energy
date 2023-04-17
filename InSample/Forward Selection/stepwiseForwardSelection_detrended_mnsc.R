# This script creates following outputs : 
# 15 files (8 LHS, 2 prediction manner (4wk, 8wk)) except for DOilVol_4wk). 
# Each file saves the regression coefficients, t-stats, P-values, etc. 
# The orders of the variables from top to bottom in the words files
# indicate the orders of being selected in the corresponding forward selection process.


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
library(imputeTS)

# source("/if/appl/R/Functions/IFfunctions.r") # loads necessary packages and IF helper functions

# Set working directory
dir <- "file_path/InSample/Forward Selection/"
setwd(dir)

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("file_path/data/transformed_data_prices_v19.dta") 

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
                                              sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
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
data.physical <- read.dta13("file_path/data/transformed_data_physical_v19.dta") 
# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t4, DInv_t8, DProd_t4, DProd_t8, FutRet, DSpot, DOilVol, xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq,  PCAall, BEME, Mom, BasMom, 
                                                  DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "ovx_diff", 
           "vix_diff", "sdf_fullSample", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Forward selection function                                                                                     
# ------------------------------------------------------------------------

fwdSelection <- function(subset, depVar, horizon){
  
  results <- list()
  
  # Regress all variables on lagged dependent variables 
  detrended.lagDepVar.adjusted <- subset
  for (i in 1:ncol(subset)) {
    subset[,i] <- na_locf(subset[,i], na_remaining = "keep") # Impute missing values with previous week's value. (This is necessary for the purpose of dealing with stanbaugh bias in Monte_Carlo_forwardSelection. We do the same here consistency.
    
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
  final.df.less.risk <- dplyr::select(final.df, -c("ovx_diff", "sdf_fullSample"))
  
  subset <- as.matrix(final.df.less.risk)
  
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
  
  # Risk premia regressions
  #  VIX <- dplyr::select(final.df,c(VIX))
  #  risk.df <- merge(lm.df, VIX, by=0, all.x = T)
  risk.df <- dplyr::select(final.df, c(c(depVar), colnames(lm.df)[-1], c("VIX")))
  #select(final.df, c(c(depVar), colnames(lm.df)[-1], c("VIX")))
  colnames(risk.df) <- c("y", colnames(risk.df)[-1])
  results$lm.1 <- lm(y ~ ., data = risk.df)
  
  risk.df <- dplyr::select(final.df, c(c(depVar), colnames(lm.df)[-1], c("vix_diff")))
  colnames(risk.df) <- c("y", colnames(risk.df)[-1])
  results$lm.2 <- lm(y ~ ., data = risk.df)
  
  risk.df <- dplyr::select(final.df, c(c(depVar), colnames(lm.df)[-1], c("ovx_diff")))
  colnames(risk.df) <- c("y", colnames(risk.df)[-1])
  results$lm.3 <- lm(y ~ ., data = risk.df)
  
  risk.df <- dplyr::select(final.df, c(c(depVar), colnames(lm.df)[-1], c("sdf_fullSample")))
  colnames(risk.df) <- c("y", colnames(risk.df)[-1])
  results$lm.4 <- lm(y ~ ., data = risk.df)
  
  return(results)
}

# ------------------------------------------------------------------------
# Futures Returns - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t4, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
FutRet_t4.0 <- results$lm.0
FutRet_t4.1 <- results$lm.1
FutRet_t4.2 <- results$lm.2
FutRet_t4.3 <- results$lm.3
FutRet_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "FutRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
FutRet_t8.0 <- results$lm.0
FutRet_t8.1 <- results$lm.1
FutRet_t8.2 <- results$lm.2
FutRet_t8.3 <- results$lm.3
FutRet_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
# Spot Price Change - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t4, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_4wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
DSpot_t4.0 <- results$lm.0
DSpot_t4.1 <- results$lm.1
DSpot_t4.2 <- results$lm.2
DSpot_t4.3 <- results$lm.3
DSpot_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DSpot_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DSpot_t8.0 <- results$lm.0
DSpot_t8.1 <- results$lm.1
DSpot_t8.2 <- results$lm.2
DSpot_t8.3 <- results$lm.3
DSpot_t8.4 <- results$lm.4

# We don't include DOilVol regression in 4-week horizon because Bloomberg evaluates it in 30-day windows 
# and thus overlaps with predictors in terms of timing
# ------------------------------------------------------------------------
# Oil Volatility - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
# subset <- dplyr::select(detrended.prices,
#                        c(DOilVol_t4, 
#                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
#                          sp500Ret, StkIdx, basis, WIPI_4wk, artcount, entropy, sCo, fCo, 
#                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
#                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
#                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

# depVar <- "DOilVol_t4"
# horizon <- 4

# results <- fwdSelection(subset, depVar, horizon)
# DOilVol_t4.0 <- results$lm.0
# DOilVol_t4.1 <- results$lm.1
# DOilVol_t4.2 <- results$lm.2
# DOilVol_t4.3 <- results$lm.3
# DOilVol_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DOilVol_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DOilVol_t8.0 <- results$lm.0
DOilVol_t8.1 <- results$lm.1
DOilVol_t8.2 <- results$lm.2
DOilVol_t8.3 <- results$lm.3
DOilVol_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
# Exxon stock returns - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(xomRet_t4, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
xomRet_t4.0 <- results$lm.0
xomRet_t4.1 <- results$lm.1
xomRet_t4.2 <- results$lm.2
xomRet_t4.3 <- results$lm.3
xomRet_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "xomRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
xomRet_t8.0 <- results$lm.0
xomRet_t8.1 <- results$lm.1
xomRet_t8.2 <- results$lm.2
xomRet_t8.3 <- results$lm.3
xomRet_t8.4 <- results$lm.4

# ------------------------------------------------------------------------
# BP stock returns - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t4, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
bpRet_t4.0 <- results$lm.0
bpRet_t4.1 <- results$lm.1
bpRet_t4.2 <- results$lm.2
bpRet_t4.3 <- results$lm.3
bpRet_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "bpRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
bpRet_t8.0 <- results$lm.0
bpRet_t8.1 <- results$lm.1
bpRet_t8.2 <- results$lm.2
bpRet_t8.3 <- results$lm.3
bpRet_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
# Shell stock returns - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t4, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
rdsaRet_t4.0 <- results$lm.0
rdsaRet_t4.1 <- results$lm.1
rdsaRet_t4.2 <- results$lm.2
rdsaRet_t4.3 <- results$lm.3
rdsaRet_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "rdsaRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
rdsaRet_t8.0 <- results$lm.0
rdsaRet_t8.1 <- results$lm.1
rdsaRet_t8.2 <- results$lm.2
rdsaRet_t8.3 <- results$lm.3
rdsaRet_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
# Oil inventories - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t4, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
DInv_t4.0 <- results$lm.0
DInv_t4.1 <- results$lm.1
DInv_t4.2 <- results$lm.2
DInv_t4.3 <- results$lm.3
DInv_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DInv_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DInv_t8.0 <- results$lm.0
DInv_t8.1 <- results$lm.1
DInv_t8.2 <- results$lm.2
DInv_t8.3 <- results$lm.3
DInv_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
# Oil production - 4-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DProd_t4, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t4"
horizon <- 4

results <- fwdSelection(subset, depVar, horizon)
DProd_t4.0 <- results$lm.0
DProd_t4.1 <- results$lm.1
DProd_t4.2 <- results$lm.2
DProd_t4.3 <- results$lm.3
DProd_t4.4 <- results$lm.4


# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

depVar <- "DProd_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DProd_t8.0 <- results$lm.0
DProd_t8.1 <- results$lm.1
DProd_t8.2 <- results$lm.2
DProd_t8.3 <- results$lm.3
DProd_t8.4 <- results$lm.4


# ------------------------------------------------------------------------
#  Output regression results using stepwise forward selection                                                                              
# ------------------------------------------------------------------------

# Futures returns forward selection regressions results with staggered risk premia measures
tab_model(FutRet_t8.0, FutRet_t8.1, FutRet_t8.2, FutRet_t8.3, FutRet_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01), title = "Predicting 8-weeks ahead Oil Futures Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_FutRet_8wks.doc")

# Spot price changes forward selection regressions results with staggered risk premia measures
tab_model(DSpot_t8.0, DSpot_t8.1, DSpot_t8.2, DSpot_t8.3, DSpot_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01), title = "Predicting 8-weeks ahead Spot Price Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DSpot_8wks.doc")

# Oil price volatility changes forward selection regressions results with staggered risk premia measures
tab_model(DOilVol_t8.0, DOilVol_t8.1, DOilVol_t8.2, DOilVol_t8.3, DOilVol_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Volatility, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DOilVol_8wks.doc")

# Exxon stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(xomRet_t8.0, xomRet_t8.1, xomRet_t8.2, xomRet_t8.3, xomRet_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Exxon Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_xomRet_8wks.doc")

# BP stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(bpRet_t8.0, bpRet_t8.1, bpRet_t8.2, bpRet_t8.3, bpRet_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead BP Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_bpRet_8wks.doc")

# Royal Dutch Shell stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(rdsaRet_t8.0, rdsaRet_t8.1, rdsaRet_t8.2, rdsaRet_t8.3, rdsaRet_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Royal Dutch Shell Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_rdsaRet_8wks.doc")

# Oil inventory changes forward selection regressions results with staggered risk premia measures
tab_model(DInv_t8.0, DInv_t8.1, DInv_t8.2, DInv_t8.3, DInv_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Inventory Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DInv_8wks.doc")

# Oil production changes forward selection regressions results with staggered risk premia measures
tab_model(DProd_t8.0, DProd_t8.1, DProd_t8.2, DProd_t8.3, DProd_t8.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Production Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DProd_8wks.doc")


############## 4 Week Horizon #################
# Futures returns forward selection regressions results with staggered risk premia measures
tab_model(FutRet_t4.0, FutRet_t4.1, FutRet_t4.2, FutRet_t4.3, FutRet_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Oil Futures Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_FutRet_4wks.doc")

# Spot price changes forward selection regressions results with staggered risk premia measures
tab_model(DSpot_t4.0, DSpot_t4.1, DSpot_t4.2, DSpot_t4.3, DSpot_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Spot Price Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DSpot_4wks.doc")

# Oil price volatility changes forward selection regressions results with staggered risk premia measures
# tab_model(DOilVol_t4.0, DOilVol_t4.1, DOilVol_t4.2, DOilVol_t4.3, DOilVol_t4.4,
#           show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Oil Volatility, including Risk Premia Measures",
#           vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DOilVol_4wks.doc")

# Exxon stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(xomRet_t4.0, xomRet_t4.1, xomRet_t4.2, xomRet_t4.3, xomRet_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Exxon Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_xomRet_4wks.doc")

# BP stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(bpRet_t4.0, bpRet_t4.1, bpRet_t4.2, bpRet_t4.3, bpRet_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead BP Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_bpRet_4wks.doc")

# Royal Dutch Shell stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(rdsaRet_t4.0, rdsaRet_t4.1, rdsaRet_t4.2, rdsaRet_t4.3, rdsaRet_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Royal Dutch Shell Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_rdsaRet_4wks.doc")

# Oil inventory changes forward selection regressions results with staggered risk premia measures
tab_model(DInv_t4.0, DInv_t4.1, DInv_t4.2, DInv_t4.3, DInv_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Oil Inventory Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DInv_4wks.doc")

# Oil production changes forward selection regressions results with staggered risk premia measures
tab_model(DProd_t4.0, DProd_t4.1, DProd_t4.2, DProd_t4.3, DProd_t4.4,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 4-weeks ahead Oil Production Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "mnsc/forwardSelection_DProd_4wks.doc")

