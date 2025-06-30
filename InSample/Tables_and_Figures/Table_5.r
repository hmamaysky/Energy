# This script creates following outputs : 
# 8 files (8 LHS, 1 prediction manner (8wk)). 
# Each file saves the regression coefficients, t-stats, P-values, etc. 
# The base regression (first column) is produced by taking forward selection on non-text variables.
# For each subsequent column, we add one text variable at a time to the base regression.
# We want to compare the R-squared values in the end.

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
dir <- "C:/Users/13917/Desktop/CBS/Energy-RA/Table 5"
setwd(dir)

# ------------------------------------------------------------------------
# Create detrended, 4- and 8-week lagged data frames (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("C:/Users/13917/Desktop/CBS/Energy-RA/Fwd_ In sample results/transformed_data_prices_v19.3_mod_restricted_sample.dta")

# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

# Add all variables
subset.prices <- dplyr::select(data.prices, c(FutRet_t4, FutRet_t8, DSpot_t4, DSpot_t8, DOilVol_t4, DOilVol_t8, 
                                              xomRet_t4, xomRet_t8, bpRet_t4, bpRet_t8, rdsaRet_t4, rdsaRet_t8,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_4wk, WIPI_8wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                                              BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

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
data.physical <- read.dta13("C:/Users/13917/Desktop/CBS/Energy-RA/Fwd_ In sample results/transformed_data_physical_v19.3_mod_restricted_sample.dta")
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
                                                  DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

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
           "vix_diff", "sdf_fullSample", "PCAsent", "PCAfreq", "PCAall","BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt", "tone", "opec", "economic_growth", "product_prices", "gov_budget", "oil_drilling", "oil", "opu_index", "tosi")

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
  
  all_text_varlist <- c("artcount", "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", "fEpg", "sBbl", "fBbl", 
  "sRpc", "fRpc", "sEp", "fEp", "PCAsent", "PCAfreq", "PCAall", 
  "tone", "opec", "economic_growth", "product_prices", "gov_budget", "oil_drilling", "oil", "opu_index", "tosi") # External text variables

  # Remove text variables and risk premia measures from final.df

  final.df.less.text <- dplyr::select(final.df, -all_of(all_text_varlist), -c("ovx_diff", "sdf_fullSample"))

  subset <- as.matrix(final.df.less.text)

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

  # Text variables regressions
  # Add text variables one at a time to the regression
  index <- 1
  for (text_var in all_text_varlist) {
    text.df <- dplyr::select(final.df, c(c(depVar), colnames(lm.df)[-1], text_var))
    colnames(text.df) <- c("y", colnames(text.df)[-1])
    results[[ paste0("lm.", index) ]] <- lm(y ~ ., data = text.df)
    index <- index + 1
  }
  
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
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "FutRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
FutRet_t8.0 <- results$lm.0
FutRet_t8.1 <- results$lm.1
FutRet_t8.2 <- results$lm.2
FutRet_t8.3 <- results$lm.3
FutRet_t8.4 <- results$lm.4
FutRet_t8.5 <- results$lm.5
FutRet_t8.6 <- results$lm.6
FutRet_t8.7 <- results$lm.7
FutRet_t8.8 <- results$lm.8
FutRet_t8.9 <- results$lm.9
FutRet_t8.10 <- results$lm.10
FutRet_t8.11 <- results$lm.11
FutRet_t8.12 <- results$lm.12
FutRet_t8.13 <- results$lm.13
FutRet_t8.14 <- results$lm.14
FutRet_t8.15 <- results$lm.15
FutRet_t8.16 <- results$lm.16
FutRet_t8.17 <- results$lm.17
FutRet_t8.18 <- results$lm.18
FutRet_t8.19 <- results$lm.19
FutRet_t8.20 <- results$lm.20
FutRet_t8.21 <- results$lm.21
FutRet_t8.22 <- results$lm.22
FutRet_t8.23 <- results$lm.23
FutRet_t8.24 <- results$lm.24
FutRet_t8.25 <- results$lm.25
FutRet_t8.26 <- results$lm.26
FutRet_t8.27 <- results$lm.27
FutRet_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "DSpot_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DSpot_t8.0 <- results$lm.0
DSpot_t8.1 <- results$lm.1
DSpot_t8.2 <- results$lm.2
DSpot_t8.3 <- results$lm.3
DSpot_t8.4 <- results$lm.4
DSpot_t8.5 <- results$lm.5
DSpot_t8.6 <- results$lm.6
DSpot_t8.7 <- results$lm.7
DSpot_t8.8 <- results$lm.8
DSpot_t8.9 <- results$lm.9
DSpot_t8.10 <- results$lm.10
DSpot_t8.11 <- results$lm.11
DSpot_t8.12 <- results$lm.12
DSpot_t8.13 <- results$lm.13
DSpot_t8.14 <- results$lm.14
DSpot_t8.15 <- results$lm.15
DSpot_t8.16 <- results$lm.16
DSpot_t8.17 <- results$lm.17
DSpot_t8.18 <- results$lm.18
DSpot_t8.19 <- results$lm.19
DSpot_t8.20 <- results$lm.20
DSpot_t8.21 <- results$lm.21
DSpot_t8.22 <- results$lm.22
DSpot_t8.23 <- results$lm.23
DSpot_t8.24 <- results$lm.24
DSpot_t8.25 <- results$lm.25
DSpot_t8.26 <- results$lm.26
DSpot_t8.27 <- results$lm.27
DSpot_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(DOilVol_t8, 
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample,
                          PCAsent, PCAfreq, PCAall,BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "DOilVol_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DOilVol_t8.0 <- results$lm.0
DOilVol_t8.1 <- results$lm.1
DOilVol_t8.2 <- results$lm.2
DOilVol_t8.3 <- results$lm.3
DOilVol_t8.4 <- results$lm.4
DOilVol_t8.5 <- results$lm.5
DOilVol_t8.6 <- results$lm.6
DOilVol_t8.7 <- results$lm.7
DOilVol_t8.8 <- results$lm.8
DOilVol_t8.9 <- results$lm.9
DOilVol_t8.10 <- results$lm.10
DOilVol_t8.11 <- results$lm.11
DOilVol_t8.12 <- results$lm.12
DOilVol_t8.13 <- results$lm.13
DOilVol_t8.14 <- results$lm.14
DOilVol_t8.15 <- results$lm.15
DOilVol_t8.16 <- results$lm.16
DOilVol_t8.17 <- results$lm.17
DOilVol_t8.18 <- results$lm.18
DOilVol_t8.19 <- results$lm.19
DOilVol_t8.20 <- results$lm.20
DOilVol_t8.21 <- results$lm.21
DOilVol_t8.22 <- results$lm.22
DOilVol_t8.23 <- results$lm.23
DOilVol_t8.24 <- results$lm.24
DOilVol_t8.25 <- results$lm.25
DOilVol_t8.26 <- results$lm.26
DOilVol_t8.27 <- results$lm.27
DOilVol_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Exxon stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t8, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "xomRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
xomRet_t8.0 <- results$lm.0
xomRet_t8.1 <- results$lm.1
xomRet_t8.2 <- results$lm.2
xomRet_t8.3 <- results$lm.3
xomRet_t8.4 <- results$lm.4
xomRet_t8.5 <- results$lm.5
xomRet_t8.6 <- results$lm.6
xomRet_t8.7 <- results$lm.7
xomRet_t8.8 <- results$lm.8
xomRet_t8.9 <- results$lm.9
xomRet_t8.10 <- results$lm.10
xomRet_t8.11 <- results$lm.11
xomRet_t8.12 <- results$lm.12
xomRet_t8.13 <- results$lm.13
xomRet_t8.14 <- results$lm.14
xomRet_t8.15 <- results$lm.15
xomRet_t8.16 <- results$lm.16
xomRet_t8.17 <- results$lm.17
xomRet_t8.18 <- results$lm.18
xomRet_t8.19 <- results$lm.19
xomRet_t8.20 <- results$lm.20
xomRet_t8.21 <- results$lm.21
xomRet_t8.22 <- results$lm.22
xomRet_t8.23 <- results$lm.23
xomRet_t8.24 <- results$lm.24
xomRet_t8.25 <- results$lm.25
xomRet_t8.26 <- results$lm.26
xomRet_t8.27 <- results$lm.27
xomRet_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# BP stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t8, 
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "bpRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
bpRet_t8.0 <- results$lm.0
bpRet_t8.1 <- results$lm.1
bpRet_t8.2 <- results$lm.2
bpRet_t8.3 <- results$lm.3
bpRet_t8.4 <- results$lm.4
bpRet_t8.5 <- results$lm.5
bpRet_t8.6 <- results$lm.6
bpRet_t8.7 <- results$lm.7
bpRet_t8.8 <- results$lm.8
bpRet_t8.9 <- results$lm.9
bpRet_t8.10 <- results$lm.10
bpRet_t8.11 <- results$lm.11
bpRet_t8.12 <- results$lm.12
bpRet_t8.13 <- results$lm.13
bpRet_t8.14 <- results$lm.14
bpRet_t8.15 <- results$lm.15
bpRet_t8.16 <- results$lm.16
bpRet_t8.17 <- results$lm.17
bpRet_t8.18 <- results$lm.18
bpRet_t8.19 <- results$lm.19
bpRet_t8.20 <- results$lm.20
bpRet_t8.21 <- results$lm.21
bpRet_t8.22 <- results$lm.22
bpRet_t8.23 <- results$lm.23
bpRet_t8.24 <- results$lm.24
bpRet_t8.25 <- results$lm.25
bpRet_t8.26 <- results$lm.26
bpRet_t8.27 <- results$lm.27
bpRet_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Shell stock returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(rdsaRet_t8, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "rdsaRet_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
rdsaRet_t8.0 <- results$lm.0
rdsaRet_t8.1 <- results$lm.1
rdsaRet_t8.2 <- results$lm.2
rdsaRet_t8.3 <- results$lm.3
rdsaRet_t8.4 <- results$lm.4
rdsaRet_t8.5 <- results$lm.5
rdsaRet_t8.6 <- results$lm.6
rdsaRet_t8.7 <- results$lm.7
rdsaRet_t8.8 <- results$lm.8
rdsaRet_t8.9 <- results$lm.9
rdsaRet_t8.10 <- results$lm.10
rdsaRet_t8.11 <- results$lm.11
rdsaRet_t8.12 <- results$lm.12
rdsaRet_t8.13 <- results$lm.13
rdsaRet_t8.14 <- results$lm.14
rdsaRet_t8.15 <- results$lm.15
rdsaRet_t8.16 <- results$lm.16
rdsaRet_t8.17 <- results$lm.17
rdsaRet_t8.18 <- results$lm.18
rdsaRet_t8.19 <- results$lm.19
rdsaRet_t8.20 <- results$lm.20
rdsaRet_t8.21 <- results$lm.21
rdsaRet_t8.22 <- results$lm.22
rdsaRet_t8.23 <- results$lm.23
rdsaRet_t8.24 <- results$lm.24
rdsaRet_t8.25 <- results$lm.25
rdsaRet_t8.26 <- results$lm.26
rdsaRet_t8.27 <- results$lm.27
rdsaRet_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Oil inventories - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t8, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "DInv_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DInv_t8.0 <- results$lm.0
DInv_t8.1 <- results$lm.1
DInv_t8.2 <- results$lm.2
DInv_t8.3 <- results$lm.3
DInv_t8.4 <- results$lm.4
DInv_t8.5 <- results$lm.5
DInv_t8.6 <- results$lm.6
DInv_t8.7 <- results$lm.7
DInv_t8.8 <- results$lm.8
DInv_t8.9 <- results$lm.9
DInv_t8.10 <- results$lm.10
DInv_t8.11 <- results$lm.11
DInv_t8.12 <- results$lm.12
DInv_t8.13 <- results$lm.13
DInv_t8.14 <- results$lm.14
DInv_t8.15 <- results$lm.15
DInv_t8.16 <- results$lm.16
DInv_t8.17 <- results$lm.17
DInv_t8.18 <- results$lm.18
DInv_t8.19 <- results$lm.19
DInv_t8.20 <- results$lm.20
DInv_t8.21 <- results$lm.21
DInv_t8.22 <- results$lm.22
DInv_t8.23 <- results$lm.23
DInv_t8.24 <- results$lm.24
DInv_t8.25 <- results$lm.25
DInv_t8.26 <- results$lm.26
DInv_t8.27 <- results$lm.27
DInv_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
# Oil production - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical, 
                        c(DProd_t8, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, ovx_diff, vix_diff, sdf_fullSample, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

depVar <- "DProd_t8"
horizon <- 8

results <- fwdSelection(subset, depVar, horizon)
DProd_t8.0 <- results$lm.0
DProd_t8.1 <- results$lm.1
DProd_t8.2 <- results$lm.2
DProd_t8.3 <- results$lm.3
DProd_t8.4 <- results$lm.4
DProd_t8.5 <- results$lm.5
DProd_t8.6 <- results$lm.6
DProd_t8.7 <- results$lm.7
DProd_t8.8 <- results$lm.8
DProd_t8.9 <- results$lm.9
DProd_t8.10 <- results$lm.10
DProd_t8.11 <- results$lm.11
DProd_t8.12 <- results$lm.12
DProd_t8.13 <- results$lm.13
DProd_t8.14 <- results$lm.14
DProd_t8.15 <- results$lm.15
DProd_t8.16 <- results$lm.16
DProd_t8.17 <- results$lm.17
DProd_t8.18 <- results$lm.18
DProd_t8.19 <- results$lm.19
DProd_t8.20 <- results$lm.20
DProd_t8.21 <- results$lm.21
DProd_t8.22 <- results$lm.22
DProd_t8.23 <- results$lm.23
DProd_t8.24 <- results$lm.24
DProd_t8.25 <- results$lm.25
DProd_t8.26 <- results$lm.26
DProd_t8.27 <- results$lm.27
DProd_t8.28 <- results$lm.28

# ------------------------------------------------------------------------
#  Output regression results using stepwise forward selection                                                                              
# ------------------------------------------------------------------------

# Futures returns forward selection regressions results with staggered risk premia measures
tab_model(FutRet_t8.0, FutRet_t8.1, FutRet_t8.2, FutRet_t8.3, FutRet_t8.4, FutRet_t8.5, FutRet_t8.6, FutRet_t8.7, FutRet_t8.8,
          FutRet_t8.9, FutRet_t8.10, FutRet_t8.11, FutRet_t8.12, FutRet_t8.13, FutRet_t8.14, FutRet_t8.15, FutRet_t8.16,
          FutRet_t8.17, FutRet_t8.18, FutRet_t8.19, FutRet_t8.20, FutRet_t8.21, FutRet_t8.22, FutRet_t8.23, FutRet_t8.24,
          FutRet_t8.25, FutRet_t8.26, FutRet_t8.27, FutRet_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01), title = "Predicting 8-weeks ahead Oil Futures Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_FutRet_8wks.doc")

# Spot price changes forward selection regressions results with staggered risk premia measures
tab_model(DSpot_t8.0, DSpot_t8.1, DSpot_t8.2, DSpot_t8.3, DSpot_t8.4, DSpot_t8.5, DSpot_t8.6, DSpot_t8.7, DSpot_t8.8,
          DSpot_t8.9, DSpot_t8.10, DSpot_t8.11, DSpot_t8.12, DSpot_t8.13, DSpot_t8.14, DSpot_t8.15, DSpot_t8.16,
          DSpot_t8.17, DSpot_t8.18, DSpot_t8.19, DSpot_t8.20, DSpot_t8.21, DSpot_t8.22, DSpot_t8.23, DSpot_t8.24,
          DSpot_t8.25, DSpot_t8.26, DSpot_t8.27, DSpot_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01), title = "Predicting 8-weeks ahead Spot Price Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_DSpot_8wks.doc")

# Oil price volatility changes forward selection regressions results with staggered risk premia measures
tab_model(DOilVol_t8.0, DOilVol_t8.1, DOilVol_t8.2, DOilVol_t8.3, DOilVol_t8.4, DOilVol_t8.5, DOilVol_t8.6, DOilVol_t8.7, DOilVol_t8.8,
          DOilVol_t8.9, DOilVol_t8.10, DOilVol_t8.11, DOilVol_t8.12, DOilVol_t8.13, DOilVol_t8.14, DOilVol_t8.15, DOilVol_t8.16,
          DOilVol_t8.17, DOilVol_t8.18, DOilVol_t8.19, DOilVol_t8.20, DOilVol_t8.21, DOilVol_t8.22, DOilVol_t8.23, DOilVol_t8.24,
          DOilVol_t8.25, DOilVol_t8.26, DOilVol_t8.27, DOilVol_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Volatility, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_DOilVol_8wks.doc")

# Exxon stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(xomRet_t8.0, xomRet_t8.1, xomRet_t8.2, xomRet_t8.3, xomRet_t8.4, xomRet_t8.5, xomRet_t8.6, xomRet_t8.7, xomRet_t8.8,
          xomRet_t8.9, xomRet_t8.10, xomRet_t8.11, xomRet_t8.12, xomRet_t8.13, xomRet_t8.14, xomRet_t8.15, xomRet_t8.16,
          xomRet_t8.17, xomRet_t8.18, xomRet_t8.19, xomRet_t8.20, xomRet_t8.21, xomRet_t8.22, xomRet_t8.23, xomRet_t8.24,
          xomRet_t8.25, xomRet_t8.26, xomRet_t8.27, xomRet_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Exxon Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_xomRet_8wks.doc")

# BP stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(bpRet_t8.0, bpRet_t8.1, bpRet_t8.2, bpRet_t8.3, bpRet_t8.4, bpRet_t8.5, bpRet_t8.6, bpRet_t8.7, bpRet_t8.8,
          bpRet_t8.9, bpRet_t8.10, bpRet_t8.11, bpRet_t8.12, bpRet_t8.13, bpRet_t8.14, bpRet_t8.15, bpRet_t8.16,
          bpRet_t8.17, bpRet_t8.18, bpRet_t8.19, bpRet_t8.20, bpRet_t8.21, bpRet_t8.22, bpRet_t8.23, bpRet_t8.24,
          bpRet_t8.25, bpRet_t8.26, bpRet_t8.27, bpRet_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead BP Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_bpRet_8wks.doc")

# Royal Dutch Shell stock returns changes forward selection regressions results with staggered risk premia measures
tab_model(rdsaRet_t8.0, rdsaRet_t8.1, rdsaRet_t8.2, rdsaRet_t8.3, rdsaRet_t8.4, rdsaRet_t8.5, rdsaRet_t8.6, rdsaRet_t8.7, rdsaRet_t8.8,
          rdsaRet_t8.9, rdsaRet_t8.10, rdsaRet_t8.11, rdsaRet_t8.12, rdsaRet_t8.13, rdsaRet_t8.14, rdsaRet_t8.15, rdsaRet_t8.16,
          rdsaRet_t8.17, rdsaRet_t8.18, rdsaRet_t8.19, rdsaRet_t8.20, rdsaRet_t8.21, rdsaRet_t8.22, rdsaRet_t8.23, rdsaRet_t8.24,
          rdsaRet_t8.25, rdsaRet_t8.26, rdsaRet_t8.27, rdsaRet_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Royal Dutch Shell Stock Returns, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_rdsaRet_8wks.doc")

# Oil inventory changes forward selection regressions results with staggered risk premia measures
tab_model(DInv_t8.0, DInv_t8.1, DInv_t8.2, DInv_t8.3, DInv_t8.4, DInv_t8.5, DInv_t8.6, DInv_t8.7, DInv_t8.8,
          DInv_t8.9, DInv_t8.10, DInv_t8.11, DInv_t8.12, DInv_t8.13, DInv_t8.14, DInv_t8.15, DInv_t8.16,
          DInv_t8.17, DInv_t8.18, DInv_t8.19, DInv_t8.20, DInv_t8.21, DInv_t8.22, DInv_t8.23, DInv_t8.24,
          DInv_t8.25, DInv_t8.26, DInv_t8.27, DInv_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Inventory Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_DInv_8wks.doc")

# Oil production changes forward selection regressions results with staggered risk premia measures
tab_model(DProd_t8.0, DProd_t8.1, DProd_t8.2, DProd_t8.3, DProd_t8.4, DProd_t8.5, DProd_t8.6, DProd_t8.7, DProd_t8.8,
          DProd_t8.9, DProd_t8.10, DProd_t8.11, DProd_t8.12, DProd_t8.13, DProd_t8.14, DProd_t8.15, DProd_t8.16,
          DProd_t8.17, DProd_t8.18, DProd_t8.19, DProd_t8.20, DProd_t8.21, DProd_t8.22, DProd_t8.23, DProd_t8.24,
          DProd_t8.25, DProd_t8.26, DProd_t8.27, DProd_t8.28,
          show.ci = F, show.p = F, p.style = "stars", p.threshold = c(0.1, 0.05, 0.01),  title = "Predicting 8-weeks ahead Oil Production Changes, including Risk Premia Measures",
          vcov.fun = "vcovHAC", file = "forwardSelection_DProd_8wks.doc")