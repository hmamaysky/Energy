# About this script:  
# All variables are candidates: baseline, text, PCs, and risk premia
# Per variable, we are using residuals after detrending, and then again 
# after taking out the lagged dependent variable that usually figures on
# the RHS of each regression (and so not per se the lagged dependent variable)
# Choose 7 variables using forward selection instead 
# Outputs: 
# Fig3.pdf: Figure 3 in the paper. Selected simulated distribution of R2 for 8wk Oil Futures Return and Volatility.
# Fig4.pdf: Figure 4 in the paper. Selected simulated distribution of t-stats for forward selected variables
#           with LHS as 8wk Oil Futures Return or Volatility.
# Word Files (4 in total)
# - 2 files for aggregated forward selection regression results (each for 4 or 8 week prediction period); 
# - Other 2 files for the AR(1) regression for error term distribution in Monte Carlo simulation (See online appendix for details)
# PDF Files (2 in total) : All simulated R2 and t-stats distributions as in Fig3 and Fig4.

# Clear Environment 
rm(list = ls())

# Libraries 
library(base)
library(dplyr)
library(sjPlot)
library(ggpubr)
library(Matrix)
library(pracma)
library(viridis)
library(cowplot)
library(sandwich)
library(readstata13)
library(selectiveInference)
library(reshape2)
library(ggplot2)
library(data.table)
library(imputeTS)

# source("/if/appl/R/Functions/IFfunctions.r") # loads necessary packages and IF helper functions

# Set working directory
dir <- "C:/Users/13917/Desktop/CBS/Energy-RA/Longer Horizon/longerHorizonAllExternal"
setwd(dir)

# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("C:/Users/13917/Desktop/CBS/Energy-RA/Longer Horizon/transformed_data_prices_v19.4_mod_restricted_sample.dta")

# Strip out the suffix we added for indicating when it's calculated
for ( col in 6:ncol(data.prices)){
  colnames(data.prices)[col] <-  sub("_M.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_T.*", "", colnames(data.prices)[col]) #works for Tue/Thu/ThFr
  colnames(data.prices)[col] <-  sub("_W.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_F.*", "", colnames(data.prices)[col])
  colnames(data.prices)[col] <-  sub("_m.*", "", colnames(data.prices)[col])
}

# Add the external tosi variable
subset.prices <- dplyr::select(data.prices, c(FutRet_t24, FutRet_t52, DSpot_t24, DSpot_t52, DOilVol_t24, DOilVol_t52, 
                                              xomRet_t24, xomRet_t52, bpRet_t24, bpRet_t52, rdsaRet_t24, rdsaRet_t52,
                                              FutRet, DSpot, DOilVol,  OilVol, DInv, xomRet, bpRet, rdsaRet,
                                              DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, WIPI_52wk, trend, artcount, 
                                              entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                                              sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
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
# Create detrended data frame (Tuesday EOP)     
# ------------------------------------------------------------------------
# Dataset for DInv and DProd
# Includes risk premia measures and VIX

data.physical <- read.dta13("C:/Users/13917/Desktop/CBS/Energy-RA/Longer Horizon/transformed_data_physical_v19.4_mod_restricted_sample.dta")

# Strip out the suffix we added for indicating when it's calculated
for ( col in 3:ncol(data.physical)){
  colnames(data.physical)[col] <-  sub("_T.*", "", colnames(data.physical)[col]) #works for both Tue and Thu
  colnames(data.physical)[col] <-  sub("_W.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_F.*", "", colnames(data.physical)[col])
  colnames(data.physical)[col] <-  sub("_m.*", "", colnames(data.physical)[col])
}

subset.physical <- dplyr::select(data.physical, c(DInv_t24, DInv_t52, DProd_t24, DProd_t52, FutRet, DSpot, DOilVol,xomRet, bpRet, rdsaRet,
                                                  OilVol, DInv, DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, WIPI_52wk, trend, artcount, 
                                                  entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, fRpc, sEp, fEp, 
                                                  VIX, vix_diff, PCAsent, PCAfreq, PCAall, BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

# Detrend all series 
detrended.physical <- subset.physical
for (i in 1:ncol(subset.physical)) {
  lm <- lm(subset.physical[,i] ~ trend, data = subset.physical, na.action=na.exclude)
  pred.value <- predict.lm(lm)
  detrended.physical[,i] <- subset.physical[,i] - pred.value
}

# Drop trend from data frame after detrending
detrended.physical <- dplyr::select(detrended.physical, -c(trend))

toLag <- c("FutRet", "DSpot", "DOilVol", "OilVol", 
           "DInv", "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "artcount", "xomRet", "bpRet", "rdsaRet",
           "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", 
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt", "tone", "opec", "economic_growth", "product_prices", "gov_budget", "oil_drilling", "oil", "opu_index", "tosi")

# ------------------------------------------------------------------------
# Forward selection function                                                                                     
# ------------------------------------------------------------------------
fwdSelection <- function(subset){
  
  results.sub <- list()
  
  subset <- as.matrix(subset)
  
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
  
  results.sub$lm <- lm(y ~ ., data = lm.df)
  results.sub$selected <- selected
  results.sub$sds <- c(sd(lm.df[,1]))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,2])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,3])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,4])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,5])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,6])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,7])))
  results.sub$sds <- c(results.sub$sds, c(sd(lm.df[,8])))
  return(results.sub)
}

# ------------------------------------------------------------------------
# MC function with forward selection
# ------------------------------------------------------------------------
# The lagged dependent variable need to be the second variable in the 
# subset data frame, needs improvement to automate this process

MC.FS <- function(subset, seed, horizon, title, sims){
  
  results <- list()
  
  # Regress all variables on lagged dependent variables 
  detrended.lagDepVar.adjusted <- subset
  for (i in 1:ncol(subset)) {
    subset[,i] <- na_locf(subset[,i], na_remaining = "keep") # Impute missing values with previous week's value. (This is necessary for the purpose of dealing with stanbaugh bias)
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
  
  # Output regression of dependent variable on lagged dependent variable for "AR1" coefficient
  results$depVar.lagDepVar <- lm(subset[,1] ~ dplyr::lag(subset[,2], horizon), data = subset)
  
  # Lag detrended and lagged-dependent-variable adjusted independent variables for forward selection
  final.df <- detrended.lagDepVar.adjusted
  for(i in 2:ncol(detrended.lagDepVar.adjusted)){
    var <- colnames(detrended.lagDepVar.adjusted)[i]
    if(sum(var == toLag) == 1){
      final.df[,i] <- dplyr::lag(detrended.lagDepVar.adjusted[,i], horizon)
    }
  }
  
  # Prep: estimate AR(h) model for dependent variable and get var(e)
  if(horizon == 4){
    AR <- lm(final.df[[1]] ~ 
               dplyr::lag(final.df[[1]], n = 1) + dplyr::lag(final.df[[1]], n = 2) + dplyr::lag(final.df[[1]], n = 3) + dplyr::lag(final.df[[1]], n = 4),
             data = final.df)
    
  } else if(horizon == 8){
    AR <- lm(final.df[[1]] ~ 
               dplyr::lag(final.df[[1]], n = 1) + dplyr::lag(final.df[[1]], n = 2) + dplyr::lag(final.df[[1]], n = 3) + dplyr::lag(final.df[[1]], n = 4) +
               dplyr::lag(final.df[[1]], n = 5) + dplyr::lag(final.df[[1]], n = 6) + dplyr::lag(final.df[[1]], n = 7) + dplyr::lag(final.df[[1]], n = 8),
             data = final.df)
    
  } else if(horizon == 24){
    AR <- lm(final.df[[1]] ~ 
              dplyr::lag(final.df[[1]], n = 1) + dplyr::lag(final.df[[1]], n = 2) + dplyr::lag(final.df[[1]], n = 3) + dplyr::lag(final.df[[1]], n = 4) +
              dplyr::lag(final.df[[1]], n = 5) + dplyr::lag(final.df[[1]], n = 6) + dplyr::lag(final.df[[1]], n = 7) + dplyr::lag(final.df[[1]], n = 8) +
              dplyr::lag(final.df[[1]], n = 9) + dplyr::lag(final.df[[1]], n = 10) + dplyr::lag(final.df[[1]], n = 11) + dplyr::lag(final.df[[1]], n = 12) +
              dplyr::lag(final.df[[1]], n = 13) + dplyr::lag(final.df[[1]], n = 14) + dplyr::lag(final.df[[1]], n = 15) + dplyr::lag(final.df[[1]], n = 16) +
              dplyr::lag(final.df[[1]], n = 17) + dplyr::lag(final.df[[1]], n = 18) + dplyr::lag(final.df[[1]], n = 19) + dplyr::lag(final.df[[1]], n = 20) +
              dplyr::lag(final.df[[1]], n = 21) + dplyr::lag(final.df[[1]], n = 22) + dplyr::lag(final.df[[1]], n = 23) + dplyr::lag(final.df[[1]], n = 24),
           data = final.df)
  } else if(horizon == 52){
    AR <- lm(final.df[[1]] ~
              dplyr::lag(final.df[[1]], n = 1) + dplyr::lag(final.df[[1]], n = 2) + dplyr::lag(final.df[[1]], n = 3) + dplyr::lag(final.df[[1]], n = 4) +
              dplyr::lag(final.df[[1]], n = 5) + dplyr::lag(final.df[[1]], n = 6) + dplyr::lag(final.df[[1]], n = 7) + dplyr::lag(final.df[[1]], n = 8) +
              dplyr::lag(final.df[[1]], n = 9) + dplyr::lag(final.df[[1]], n = 10) + dplyr::lag(final.df[[1]], n = 11) + dplyr::lag(final.df[[1]], n = 12) +
              dplyr::lag(final.df[[1]], n = 13) + dplyr::lag(final.df[[1]], n = 14) + dplyr::lag(final.df[[1]], n = 15) + dplyr::lag(final.df[[1]], n = 16) +
              dplyr::lag(final.df[[1]], n = 17) + dplyr::lag(final.df[[1]], n = 18) + dplyr::lag(final.df[[1]], n = 19) + dplyr::lag(final.df[[1]], n = 20) +
              dplyr::lag(final.df[[1]], n = 21) + dplyr::lag(final.df[[1]], n = 22) + dplyr::lag(final.df[[1]], n = 23) + dplyr::lag(final.df[[1]], n = 24) +
              dplyr::lag(final.df[[1]], n = 25) + dplyr::lag(final.df[[1]], n = 26) + dplyr::lag(final.df[[1]], n = 27) + dplyr::lag(final.df[[1]], n = 28) +
              dplyr::lag(final.df[[1]], n = 29) + dplyr::lag(final.df[[1]], n = 30) + dplyr::lag(final.df[[1]], n = 31) + dplyr::lag(final.df[[1]], n = 32) +
              dplyr::lag(final.df[[1]], n = 33) + dplyr::lag(final.df[[1]], n = 34) + dplyr::lag(final.df[[1]], n = 35) + dplyr::lag(final.df[[1]], n = 36) +
              dplyr::lag(final.df[[1]], n = 37) + dplyr::lag(final.df[[1]], n = 38) + dplyr::lag(final.df[[1]], n = 39) + dplyr::lag(final.df[[1]], n = 40) +
              dplyr::lag(final.df[[1]], n = 41) + dplyr::lag(final.df[[1]], n = 42) + dplyr::lag(final.df[[1]], n = 43) + dplyr::lag(final.df[[1]], n = 44) +
              dplyr::lag(final.df[[1]], n = 45) + dplyr::lag(final.df[[1]], n = 46) + dplyr::lag(final.df[[1]], n = 47) + dplyr::lag(final.df[[1]], n = 48) +
              dplyr::lag(final.df[[1]], n = 49) + dplyr::lag(final.df[[1]], n = 50) + dplyr::lag(final.df[[1]], n = 51) + dplyr::lag(final.df[[1]], n = 52),
           data = final.df)
  } else {
    print("Error! Check horizon for AR.")
  }
  
  y_hat <- predict.lm(AR, final.df)
  y_resid1 <- final.df[[1]] - y_hat
  # get corr(y_resid1 and detrended.lagDepVar.adjusted[]). we should not use final.df because the rhs are lagged for fwdSelection.
  cor_all <- array()
  for (i in 3:ncol(detrended.lagDepVar.adjusted)) {
    cor_all[i-2] <- abs(cor(as.data.frame(detrended.lagDepVar.adjusted[[i]]),as.data.frame(y_resid1), use = "pairwise.complete.obs"))
  }
  cor_all <- as.data.frame(cor_all)
  rownames(cor_all) <- colnames(detrended.lagDepVar.adjusted)[3:ncol(detrended.lagDepVar.adjusted)]
  # then get five most correlated ones
  cor_top5 <- tail(arrange(cor_all,cor_all),5)
  # extract top 5 from detrended.lagDepVar.adjusted
  top5list <- rownames(cor_top5)
  idx <- match(top5list,names(detrended.lagDepVar.adjusted))
  val_top5 <- detrended.lagDepVar.adjusted[,idx]
  new.df <- cbind(as.data.frame(y_resid1),val_top5)
  fit_new <- lm(y_resid1~., data=new.df)
  y_hat_resid1 <- predict.lm(fit_new, new.df)
  y_resid2 <- new.df[[1]] -y_hat_resid1
  
  var_e <- var(na.omit(y_resid2))
  
  # Prep for bootstrap
  fs.base <- fwdSelection(final.df)
  results$fs.reg <- fs.base$lm
  results$sd <- fs.base$sds
  
  rsq.base <- summary(fs.base$lm)$adj.r.squared
  rsq.sims <- matrix(nrow = sims, ncol = 1)
  
  tstat.base <- summary(fs.base$lm)$coefficient[2:8,3]
  tstat.sims <- matrix(nrow = sims, ncol = 7)
  tstat.vars <- matrix(nrow = sims, ncol = 7)
  tstat.coef <- matrix(nrow = sims, ncol = 7)
  
  # Run the boostrap
  set.seed(seed)
  
  for(m in 1:sims){
    # 1. Set y_1:h, ... y_K:K+h-1 = 0
    if(horizon == 4){
      if(rowSums(is.na(val_top5[5,]))>0){
        obs <- sum(!is.na(subset[1])) + 101 # It's not 104 because we lose three, i.e., three missing values in DProd for obs 5, 6, and 7. We thus use rhs variable from row 8 if horizon =4.
        val_top5 <- val_top5[-c(5:7),]
      } else {
        obs <- sum(!is.na(subset[1])) + 104
      }
      y <- c(rep(NA, obs))
      y[1:4] <- 0
    } else if(horizon == 8){
      obs <- sum(!is.na(subset[1])) + 108
      y <- c(rep(NA, obs))
      y[1:8] <- 0
    } else if(horizon == 24){
      obs <- sum(!is.na(subset[1])) + 124
      y <- c(rep(NA, obs))
      y[1:24] <- 0
    } else if(horizon == 52){
      obs <- sum(!is.na(subset[1])) + 152
      y <- c(rep(NA, obs))
      y[1:52] <- 0
    } else {
      print("Error! Check step 1 of bootstrap.")
    }
    
    # 2. Draw e_k+1:k+h from a normal distribution with mean zero and variance var(e)
    e <- rep(NA, obs)
    e <- rnorm(obs, mean = 0, sd = sqrt(var_e))
    
    # +) Create final X matrix for the use of stambaugh bias.
    x_top5 <-  cbind(const=1,val_top5)
    zeros <- data.frame(matrix(0, ncol=6, nrow=100))
    colnames(zeros) <- colnames(x_top5)
    x_top5 <- rbind(zeros,x_top5)
    # 3. Use the AR(8) relationship to generate y_k+1:k+h
    # 4. Run the model for 100 steps as a burn-in period
    # 5. 
    if(horizon == 4){
      coefs <- summary(AR)$coefficients[2:5,1]
      x_top5[1:104,] <- 0
      
      for(i in 5:obs){
        
        y[i] <- coefs[1]*y[i-1] + coefs[2]*y[i-2] + coefs[3]*y[i-3] + coefs[4]*y[i-4] + e[i] + 
          sum(x_top5[i,]*fit_new$coefficient)
      }
      
      y_sim <- y[105:obs]
      
    } else if(horizon == 8){
      coefs <- summary(AR)$coefficients[2:9]
      x_top5[1:108,] <- 0
      
      for(i in 9:obs){
        y[i] <- coefs[1]*y[i-1] + coefs[2]*y[i-2] + coefs[3]*y[i-3] + coefs[4]*y[i-4] +
          coefs[5]*y[i-5] + coefs[6]*y[i-6] + coefs[7]*y[i-7] + coefs[8]*y[i-8] + e[i] + 
          sum(x_top5[i,]*fit_new$coefficient)
      }
      
      y_sim <- y[109:obs]
      
    } else if(horizon == 24) {
      coefs <- summary(AR)$coefficients[2:25]
      x_top5[1:124, ] <- 0

      for (i in 25:obs) {
          y[i] <- coefs[1]*y[i-1]  + coefs[2]*y[i-2]  + coefs[3]*y[i-3]  + coefs[4]*y[i-4] +
                  coefs[5]*y[i-5]  + coefs[6]*y[i-6]  + coefs[7]*y[i-7]  + coefs[8]*y[i-8] +
                  coefs[9]*y[i-9]  + coefs[10]*y[i-10]+ coefs[11]*y[i-11]+ coefs[12]*y[i-12] +
                  coefs[13]*y[i-13]+ coefs[14]*y[i-14]+ coefs[15]*y[i-15]+ coefs[16]*y[i-16] +
                  coefs[17]*y[i-17]+ coefs[18]*y[i-18]+ coefs[19]*y[i-19]+ coefs[20]*y[i-20] +
                  coefs[21]*y[i-21]+ coefs[22]*y[i-22]+ coefs[23]*y[i-23]+ coefs[24]*y[i-24] +
                  e[i] + sum(x_top5[i, ] * fit_new$coefficient)
    }

      y_sim <- y[125:obs]  

    } else if (horizon == 52) {
      coefs <- summary(AR)$coefficients[2:53]
      x_top5[1:152, ] <- 0
      for (i in 53:obs) {
          y[i] <- coefs[1]*y[i-1]   + coefs[2]*y[i-2]   + coefs[3]*y[i-3]   + coefs[4]*y[i-4]  +
                  coefs[5]*y[i-5]   + coefs[6]*y[i-6]   + coefs[7]*y[i-7]   + coefs[8]*y[i-8]  +
                  coefs[9]*y[i-9]   + coefs[10]*y[i-10] + coefs[11]*y[i-11] + coefs[12]*y[i-12] +
                  coefs[13]*y[i-13] + coefs[14]*y[i-14] + coefs[15]*y[i-15] + coefs[16]*y[i-16] +
                  coefs[17]*y[i-17] + coefs[18]*y[i-18] + coefs[19]*y[i-19] + coefs[20]*y[i-20] +
                  coefs[21]*y[i-21] + coefs[22]*y[i-22] + coefs[23]*y[i-23] + coefs[24]*y[i-24] +
                  coefs[25]*y[i-25] + coefs[26]*y[i-26] + coefs[27]*y[i-27] + coefs[28]*y[i-28] +
                  coefs[29]*y[i-29] + coefs[30]*y[i-30] + coefs[31]*y[i-31] + coefs[32]*y[i-32] +
                  coefs[33]*y[i-33] + coefs[34]*y[i-34] + coefs[35]*y[i-35] + coefs[36]*y[i-36] +
                  coefs[37]*y[i-37] + coefs[38]*y[i-38] + coefs[39]*y[i-39] + coefs[40]*y[i-40] +
                  coefs[41]*y[i-41] + coefs[42]*y[i-42] + coefs[43]*y[i-43] + coefs[44]*y[i-44] +
                  coefs[45]*y[i-45] + coefs[46]*y[i-46] + coefs[47]*y[i-47] + coefs[48]*y[i-48] +
                  coefs[49]*y[i-49] + coefs[50]*y[i-50] + coefs[51]*y[i-51] + coefs[52]*y[i-52] +
                  e[i] + sum(x_top5[i, ] * fit_new$coefficient)
      }

      y_sim <- y[153:obs]

    } else {
      print("Error!")
    }
    
    # 6. 
    indep.vars <- dplyr::select(final.df, -c(1))
    # R counts observation including missing values, i.e. 1147. But we should count only nonmissing rows. ex) FutRet_t8 has 1139 nonmissing obs.
    # y_sim are built with that principle, and so we also do the same for indep.vars and match the y_sim and indep.var.
    
    indep.vars <- indep.vars %>% filter(row_number() %in% as.numeric(count(indep.vars)-length(y_sim)+1): as.numeric(count(indep.vars)))
    subset <- cbind(y_sim, indep.vars)
    fs.sims <- fwdSelection(subset)
    rsq.sims[m] <- summary(fs.sims$lm)$adj.r.squared
    tstat.sims[m,] <- summary(fs.sims$lm)$coefficients[2:8, 3]
    tstat.coef[m,] <- summary(fs.sims$lm)$coefficients[2:8, 1]
    tstat.vars[m,] <- fs.sims$selected
  }
  
  results$tstat.sims <- tstat.sims
  results$tstat.vars <- tstat.vars
  results$tstat.coef <- tstat.coef
  
  # Compute p-values 
  pval1 <- (sum(rsq.sims < rsq.base)/sims) * 100
  
  pval2 <- (sum(tstat.sims[,1] < tstat.base[1])/sims) 
  pval3 <- (sum(tstat.sims[,2] < tstat.base[2])/sims) 
  pval4 <- (sum(tstat.sims[,3] < tstat.base[3])/sims) 
  pval5 <- (sum(tstat.sims[,4] < tstat.base[4])/sims) 
  pval6 <- (sum(tstat.sims[,5] < tstat.base[5])/sims) 
  pval7 <- (sum(tstat.sims[,6] < tstat.base[6])/sims) 
  pval8 <- (sum(tstat.sims[,7] < tstat.base[7])/sims) 
  
  # Take min of p-value or (100 - p-value)
  #pval1 <- min(pval1, 100-pval1)
  pval2 <- min(pval2, 1-pval2)
  pval3 <- min(pval3, 1-pval3)
  pval4 <- min(pval4, 1-pval4)
  pval5 <- min(pval5, 1-pval5)
  pval6 <- min(pval6, 1-pval6)
  pval7 <- min(pval7, 1-pval7)
  pval8 <- min(pval8, 1-pval8)
  
  data.rsq <- as.data.table(rsq.sims)
  dens.rsq <- melt(data.rsq)
  
  data.tstat <- as.data.table(tstat.sims)
  dens.tstat <- melt(data.tstat)
  
  # Density plot for adjusted r-squared
  results$plot.rsq <- ggplot(data = dens.rsq, aes(x = value, group = variable)) + 
    geom_density(aes(fill = factor(variable, levels = c("V1"), labels = c(paste0("P-value: ",pval1,"%\nDifference: ", round((rsq.base-mean(rsq.sims))*100,2), " ppt")))), alpha = 0.3) +
    geom_vline(xintercept = mean(rsq.sims)) +
    geom_text(mapping = aes(x = mean(rsq.sims), y = 0, label = "Mean", hjust = -0.5, vjust = -1)) +
    ylab(NULL) + xlab(NULL) + 
    theme(plot.margin = unit(c(4, 6, 4, 2), "mm"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10, hjust = 0),
          plot.caption = element_text(size = 10),
          axis.ticks.length = unit(0.17, "cm"), 
          axis.line = element_line(size = 0.45),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(margin = unit(c(5,0,0,0), "pt")),
          axis.text.y = element_text(margin = unit(c(0,5,0,0), "pt")),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.width = unit(.4, "cm"), 
          legend.key.height= unit(.4, "cm"),
          legend.spacing.y = unit(1.5, "mm"),
          legend.position = c(0.85, 0.87)) +
    labs(title = paste0(title," (",horizon,"-week horizon)"), 
         subtitle = paste0("Density (N = ",sims,")"),
         caption = paste0("Note: the adjusted R-squared value in the baseline model \n is ",(round((rsq.base*100), 2)),"%."), 
         fill = "") 
  
  # Density plot for t-statistics
  results$plot.tstat <- ggplot(data = dens.tstat, aes(x = value, group = variable)) + 
    geom_density(aes(fill = factor(variable, levels = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), 
                                   labels = c(paste0("Variable 1: ",round(pval2, 3)*100), paste0("Variable 2: ",round(pval3, 3)*100), paste0("Variable 3: ",round(pval4, 3)*100), paste0("Variable 4: ",round(pval5, 3)*100),
                                              paste0("Variable 5: ",round(pval6, 3)*100), paste0("Variable 6: ",round(pval7, 3)*100), paste0("Variable 7: ",round(pval8, 3)*100)))), alpha = 0.3) +
    ylab(NULL) + xlab(NULL) + 
    theme(plot.margin = unit(c(4, 6, 4, 2), "mm"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10, hjust = 0),
          plot.caption = element_text(size = 10),
          axis.ticks.length = unit(0.17, "cm"), 
          axis.line = element_line(size = 0.45),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(margin = unit(c(5,0,0,0), "pt")),
          axis.text.y = element_text(margin = unit(c(0,5,0,0), "pt")),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.width = unit(.4, "cm"), 
          legend.key.height= unit(.4, "cm"),
          legend.spacing.y = unit(1.5, "mm"),
          legend.position = c(0.9, 0.85)) +
    labs(title = paste0(title," (",horizon,"-week horizon)"), 
         subtitle = paste0("Density (N = ",sims,")"),
         caption = paste0("Note: the t-statistics for the 7 variables chosen via forward selection are: ",
                          round(tstat.base[1], 2),", ", round(tstat.base[2], 2),", ", round(tstat.base[3], 2),", ", round(tstat.base[4], 2),", ", 
                          round(tstat.base[5], 2),", ", round(tstat.base[6], 2),", and ", round(tstat.base[7], 2),"."), 
         fill = "P-value (%)") 
  
  print(min(rsq.sims))
  print(max(rsq.sims))
  print(min(tstat.sims))
  print(max(tstat.sims))
  
  print("Mean R-Square")
  print(mean(rsq.sims))
  
  return(results)
}

# ------------------------------------------------------------------------
# Futures Returns - 24-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t24, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 876
title <- "Oil Futures Returns"

tic()
FutRet_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()

#FutRet_4wk <- FutRet_4wk + 
#  scale_x_continuous(expand = c(0,0), limits = c(0, 0.25), breaks = seq(0, 0.25, by = 0.05)) +
#  scale_y_continuous(expand = c(0,0), position = "left",
#                     limits = c(0,15), breaks = seq(0, 15, by = 5))

# ------------------------------------------------------------------------
# Futures Returns - 52-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t52, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Spot Price Change - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t24, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_24wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 2139
title <- "Spot Price Changes"

tic()
DSpot_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Spot Price Change - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t52, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_52wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                       c(DOilVol_t24, 
                         DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                         sp500Ret, StkIdx, basis, WIPI_24wk, artcount, entropy, sCo, fCo, 
                         sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                         fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                         BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 55552
title <- "Oil Volatility"

tic()
DOilVol_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Oil Volatility - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t52,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_52wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Exxon stock returns - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t24, 
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 4401223
title <- "Exxon Stock Returns"

tic()
xomRet_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Exxon stock returns - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(xomRet_t52,
                          xomRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# BP stock returns - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t24, 
                          bpRet, FutRet, DSpot, DOilVol,OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 33312220
title <- "BP Stock Returns"

tic()
bpRet_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# BP stock returns - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(bpRet_t52,
                          bpRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Shell stock returns - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(rdsaRet_t24, 
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 309912
title <- "Shell Stock Returns"

tic()
rdsaRet_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Shell stock returns - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(rdsaRet_t52,
                          rdsaRet, FutRet, DSpot, DOilVol, OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Oil inventories - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t24, 
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 751046
title <- "Oil Inventories"

tic()
DInv_4wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Oil inventories - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DInv_t52,
                          DInv, FutRet, DSpot, DOilVol,  OilVol, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Oil production - 24-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DProd_t24, 
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_24wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 24
sims <- 1000
seed <- 20227
title <- "Oil Production"

tic()
DProd_4wk<- MC.FS(subset, seed, horizon, title, sims)
toc()


# ------------------------------------------------------------------------
# Oil production - 52-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.physical,
                        c(DProd_t52,
                          DProd, FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_52wk, artcount, 
                          entropy, sCo, fCo, sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, 
                          sRpc, fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt, tone, opec, economic_growth, product_prices, gov_budget, oil_drilling, oil, opu_index, tosi))

horizon <- 52
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 52wk #####
# FutRet
FutRet_8wk$sdcoef <- FutRet_8wk$fs.reg
FutRet_8wk$sdcoef$coefficients<-c(FutRet_8wk$fs.reg$coefficients[1],FutRet_8wk$fs.reg$coefficients[2:8]/FutRet_8wk$sd[1]*FutRet_8wk$sd[-1])
# DSpot
DSpot_8wk$sdcoef <- DSpot_8wk$fs.reg
DSpot_8wk$sdcoef$coefficients<-c(DSpot_8wk$fs.reg$coefficients[1],DSpot_8wk$fs.reg$coefficients[2:8]/DSpot_8wk$sd[1]*DSpot_8wk$sd[-1])
# DOilVol
DOilVol_8wk$sdcoef <- DOilVol_8wk$fs.reg
DOilVol_8wk$sdcoef$coefficients<-c(DOilVol_8wk$fs.reg$coefficients[1],DOilVol_8wk$fs.reg$coefficients[2:8]/DOilVol_8wk$sd[1]*DOilVol_8wk$sd[-1])
# xomRet
xomRet_8wk$sdcoef <- xomRet_8wk$fs.reg
xomRet_8wk$sdcoef$coefficients<-c(xomRet_8wk$fs.reg$coefficients[1],xomRet_8wk$fs.reg$coefficients[2:8]/xomRet_8wk$sd[1]*xomRet_8wk$sd[-1])
# bpRet
bpRet_8wk$sdcoef <- bpRet_8wk$fs.reg
bpRet_8wk$sdcoef$coefficients<-c(bpRet_8wk$fs.reg$coefficients[1],bpRet_8wk$fs.reg$coefficients[2:8]/bpRet_8wk$sd[1]*bpRet_8wk$sd[-1])
# rdsaRet
rdsaRet_8wk$sdcoef <- rdsaRet_8wk$fs.reg
rdsaRet_8wk$sdcoef$coefficients<-c(rdsaRet_8wk$fs.reg$coefficients[1],rdsaRet_8wk$fs.reg$coefficients[2:8]/rdsaRet_8wk$sd[1]*rdsaRet_8wk$sd[-1])
# DInv
DInv_8wk$sdcoef <- DInv_8wk$fs.reg
DInv_8wk$sdcoef$coefficients<-c(DInv_8wk$fs.reg$coefficients[1],DInv_8wk$fs.reg$coefficients[2:8]/DInv_8wk$sd[1]*DInv_8wk$sd[-1])
# DProd
DProd_8wk$sdcoef <- DProd_8wk$fs.reg
DProd_8wk$sdcoef$coefficients<-c(DProd_8wk$fs.reg$coefficients[1],DProd_8wk$fs.reg$coefficients[2:8]/DProd_8wk$sd[1]*DProd_8wk$sd[-1])

##### 24wk #####
# FutRet
FutRet_4wk$sdcoef <- FutRet_4wk$fs.reg
FutRet_4wk$sdcoef$coefficients<-c(FutRet_4wk$fs.reg$coefficients[1],FutRet_4wk$fs.reg$coefficients[2:8]/FutRet_4wk$sd[1]*FutRet_4wk$sd[-1])
# DSpot
DSpot_4wk$sdcoef <- DSpot_4wk$fs.reg
DSpot_4wk$sdcoef$coefficients<-c(DSpot_4wk$fs.reg$coefficients[1],DSpot_4wk$fs.reg$coefficients[2:8]/DSpot_4wk$sd[1]*DSpot_4wk$sd[-1])
# DOilVol
DOilVol_4wk$sdcoef <- DOilVol_4wk$fs.reg
DOilVol_4wk$sdcoef$coefficients<-c(DOilVol_4wk$fs.reg$coefficients[1],DOilVol_4wk$fs.reg$coefficients[2:8]/DOilVol_4wk$sd[1]*DOilVol_4wk$sd[-1])
# xomRet
xomRet_4wk$sdcoef <- xomRet_4wk$fs.reg
xomRet_4wk$sdcoef$coefficients<-c(xomRet_4wk$fs.reg$coefficients[1],xomRet_4wk$fs.reg$coefficients[2:8]/xomRet_4wk$sd[1]*xomRet_4wk$sd[-1])
# bpRet
bpRet_4wk$sdcoef <- bpRet_4wk$fs.reg
bpRet_4wk$sdcoef$coefficients<-c(bpRet_4wk$fs.reg$coefficients[1],bpRet_4wk$fs.reg$coefficients[2:8]/bpRet_4wk$sd[1]*bpRet_4wk$sd[-1])
# rdsaRet
rdsaRet_4wk$sdcoef <- rdsaRet_4wk$fs.reg
rdsaRet_4wk$sdcoef$coefficients<-c(rdsaRet_4wk$fs.reg$coefficients[1],rdsaRet_4wk$fs.reg$coefficients[2:8]/rdsaRet_4wk$sd[1]*rdsaRet_4wk$sd[-1])
# DInv
DInv_4wk$sdcoef <- DInv_4wk$fs.reg
DInv_4wk$sdcoef$coefficients<-c(DInv_4wk$fs.reg$coefficients[1],DInv_4wk$fs.reg$coefficients[2:8]/DInv_4wk$sd[1]*DInv_4wk$sd[-1])
# DProd
DProd_4wk$sdcoef <- DProd_4wk$fs.reg
DProd_4wk$sdcoef$coefficients<-c(DProd_4wk$fs.reg$coefficients[1],DProd_4wk$fs.reg$coefficients[2:8]/DProd_4wk$sd[1]*DProd_4wk$sd[-1])
# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 52-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk$sdcoef, DSpot_8wk$sdcoef, DOilVol_8wk$sdcoef, xomRet_8wk$sdcoef, bpRet_8wk$sdcoef, rdsaRet_8wk$sdcoef, DInv_8wk$sdcoef, DProd_8wk$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 52-week Horizon",
          vcov.fun = "vcovHAC", file = "forwardSelection_52wks_allVars_new.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk$depVar.lagDepVar, DSpot_8wk$depVar.lagDepVar, DOilVol_8wk$depVar.lagDepVar, xomRet_8wk$depVar.lagDepVar, 
          bpRet_8wk$depVar.lagDepVar, rdsaRet_8wk$depVar.lagDepVar, DInv_8wk$depVar.lagDepVar, DProd_8wk$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "forwardSelection_52wks_AR1_new.doc")

# 24-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_4wk$sdcoef, DSpot_4wk$sdcoef, DOilVol_4wk$sdcoef, xomRet_4wk$sdcoef, bpRet_4wk$sdcoef, rdsaRet_4wk$sdcoef, DInv_4wk$sdcoef, DProd_4wk$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 24-week Horizon",
          vcov.fun = "vcovHAC", file = "forwardSelection_24wks_allVars_new.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_4wk$depVar.lagDepVar, DSpot_4wk$depVar.lagDepVar, DOilVol_4wk$depVar.lagDepVar, xomRet_4wk$depVar.lagDepVar, 
          bpRet_4wk$depVar.lagDepVar, rdsaRet_4wk$depVar.lagDepVar, DInv_4wk$depVar.lagDepVar, DProd_4wk$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "forwardSelection_24wks_AR1_new.doc")

# Adjusted R-Squared
pdf("Monte_Carlo_Exhibits_fs_rsq_new.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_4wk$plot.rsq, FutRet_8wk$plot.rsq, DSpot_4wk$plot.rsq, DSpot_8wk$plot.rsq, DOilVol_4wk$plot.rsq, DOilVol_8wk$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_4wk$plot.rsq, xomRet_8wk$plot.rsq, bpRet_4wk$plot.rsq, bpRet_8wk$plot.rsq, rdsaRet_4wk$plot.rsq, rdsaRet_8wk$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_4wk$plot.rsq, DInv_8wk$plot.rsq, DProd_4wk$plot.rsq, DProd_8wk$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("Monte_Carlo_Exhibits_fs_tstat_new.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_4wk$plot.tstat, FutRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_4wk$plot.tstat, DSpot_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_4wk$plot.tstat, DOilVol_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_4wk$plot.tstat, xomRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_4wk$plot.tstat, bpRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_4wk$plot.tstat, rdsaRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_4wk$plot.tstat, DInv_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_4wk$plot.tstat, DProd_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))
dev.off()

pdf("Fig3.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")
figure <- ggarrange(FutRet_8wk$plot.rsq, DOilVol_8wk$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("", just = "center", face = "bold", size = 14))
dev.off()

pdf("Fig4.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")
figure <- ggarrange(FutRet_8wk$plot.tstat, DOilVol_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("", just = "center", face = "bold", size = 14))
dev.off()