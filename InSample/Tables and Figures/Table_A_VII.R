# About this script:  
# This script runs monte carlo forward selection for all subperiods and produce following outputs.
# Main Output: Table A.VII "Subperiod stepwise forward selection at the eight-week horizon" in word file.
# Additional Outputs:  
# - Word files for aggregated forward selection regression results  
# - Word files for the AR(1) regression for error term distribution in Monte Carlo simulation (See online appendix for details)
# - PDF Files : All simulated R2 and t-stats distributions

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
dir <- "/Users/Economist/Dropbox/Research/ncm_research/Monte Carlo/"
setwd(dir)

# Subperiod 1 - 9


#######################
#### sub1##############
#######################

# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
# ------------------------------------------------------------------------
# Dataset for FutRet, DSpot, DOilVol, xomRet, bpRet, and rdsaRet
# Includes risk premia measures and VIX
data.prices <- read.dta13("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta") 
data.prices <- subset(data.prices, date_Fri>= as.Date("1998-04-01") & date_Fri <= as.Date("1999-11-30")) #### sub1


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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


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
  results$pval <- c(pval1, pval2, pval3, pval4, pval5, pval6, pval7, pval8)
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
  
  results$sdcoef <- results$fs.reg
  results$sdcoef$coefficients<-c(results$fs.reg$coefficients[1],results$fs.reg$coefficients[2:8]/results$sd[1]*results$sd[-1])
  
  # compile coefficients, p-val, R2, CDF information for later creating Table A.VII
  results$coeff_pval = data.table(Predictors=names(results$sdcoef$coefficients), coef=round(as.numeric(results$sdcoef$coefficients), digits=2), pval=round(results$pval,3))
  results$coeff_pval <- rbind(results$coeff_pval,list("Observations", length(results$fs.reg$model$y), ""))
  results$coeff_pval <- rbind(results$coeff_pval,list("R2 / R2 adjusted", paste0(round(summary(results$plot.rsq$plot_env$fs.base$lm)$r.squared,2), "/",
                                                                                 round(summary(results$plot.rsq$plot_env$fs.base$lm)$adj.r.squared,2)),""))
  results$coeff_pval <- rbind(results$coeff_pval,list("Mean of sim. Adj. R2", round(mean(results$plot.rsq$plot_env$rsq.sims),3), ""))
  results$coeff_pval <- rbind(results$coeff_pval,list("CDF (%)", results$plot.tstat$plot_env$pval1, ""))
  
  
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
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub1 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub1$sdcoef <- FutRet_8wk_sub1$fs.reg
FutRet_8wk_sub1$sdcoef$coefficients<-c(FutRet_8wk_sub1$fs.reg$coefficients[1],FutRet_8wk_sub1$fs.reg$coefficients[2:8]/FutRet_8wk_sub1$sd[1]*FutRet_8wk_sub1$sd[-1])
# DSpot
DSpot_8wk_sub1$sdcoef <- DSpot_8wk_sub1$fs.reg
DSpot_8wk_sub1$sdcoef$coefficients<-c(DSpot_8wk_sub1$fs.reg$coefficients[1],DSpot_8wk_sub1$fs.reg$coefficients[2:8]/DSpot_8wk_sub1$sd[1]*DSpot_8wk_sub1$sd[-1])
# DOilVol
DOilVol_8wk_sub1$sdcoef <- DOilVol_8wk_sub1$fs.reg
DOilVol_8wk_sub1$sdcoef$coefficients<-c(DOilVol_8wk_sub1$fs.reg$coefficients[1],DOilVol_8wk_sub1$fs.reg$coefficients[2:8]/DOilVol_8wk_sub1$sd[1]*DOilVol_8wk_sub1$sd[-1])
# xomRet
xomRet_8wk_sub1$sdcoef <- xomRet_8wk_sub1$fs.reg
xomRet_8wk_sub1$sdcoef$coefficients<-c(xomRet_8wk_sub1$fs.reg$coefficients[1],xomRet_8wk_sub1$fs.reg$coefficients[2:8]/xomRet_8wk_sub1$sd[1]*xomRet_8wk_sub1$sd[-1])
# bpRet
bpRet_8wk_sub1$sdcoef <- bpRet_8wk_sub1$fs.reg
bpRet_8wk_sub1$sdcoef$coefficients<-c(bpRet_8wk_sub1$fs.reg$coefficients[1],bpRet_8wk_sub1$fs.reg$coefficients[2:8]/bpRet_8wk_sub1$sd[1]*bpRet_8wk_sub1$sd[-1])
# rdsaRet
rdsaRet_8wk_sub1$sdcoef <- rdsaRet_8wk_sub1$fs.reg
rdsaRet_8wk_sub1$sdcoef$coefficients<-c(rdsaRet_8wk_sub1$fs.reg$coefficients[1],rdsaRet_8wk_sub1$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub1$sd[1]*rdsaRet_8wk_sub1$sd[-1])
# DInv
DInv_8wk_sub1$sdcoef <- DInv_8wk_sub1$fs.reg
DInv_8wk_sub1$sdcoef$coefficients<-c(DInv_8wk_sub1$fs.reg$coefficients[1],DInv_8wk_sub1$fs.reg$coefficients[2:8]/DInv_8wk_sub1$sd[1]*DInv_8wk_sub1$sd[-1])
# DProd
DProd_8wk_sub1$sdcoef <- DProd_8wk_sub1$fs.reg
DProd_8wk_sub1$sdcoef$coefficients<-c(DProd_8wk_sub1$fs.reg$coefficients[1],DProd_8wk_sub1$fs.reg$coefficients[2:8]/DProd_8wk_sub1$sd[1]*DProd_8wk_sub1$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub1$sdcoef, DSpot_8wk_sub1$sdcoef, DOilVol_8wk_sub1$sdcoef, xomRet_8wk_sub1$sdcoef, bpRet_8wk_sub1$sdcoef, rdsaRet_8wk_sub1$sdcoef, DInv_8wk_sub1$sdcoef, DProd_8wk_sub1$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wk_sub1s_allVars_new_sub1.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub1$depVar.lagDepVar, DSpot_8wk_sub1$depVar.lagDepVar, DOilVol_8wk_sub1$depVar.lagDepVar, xomRet_8wk_sub1$depVar.lagDepVar, 
          bpRet_8wk_sub1$depVar.lagDepVar, rdsaRet_8wk_sub1$depVar.lagDepVar, DInv_8wk_sub1$depVar.lagDepVar, DProd_8wk_sub1$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wk_sub1s_AR1_new_sub1.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub1.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub1$plot.rsq, DSpot_8wk_sub1$plot.rsq, DOilVol_8wk_sub1$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub1$plot.rsq, bpRet_8wk_sub1$plot.rsq, rdsaRet_8wk_sub1$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub1$plot.rsq, DProd_8wk_sub1$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub1.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub1$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()

#######################
#### sub2##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub2<- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub2 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub2$sdcoef <- FutRet_8wk_sub2$fs.reg
FutRet_8wk_sub2$sdcoef$coefficients<-c(FutRet_8wk_sub2$fs.reg$coefficients[1],FutRet_8wk_sub2$fs.reg$coefficients[2:8]/FutRet_8wk_sub2$sd[1]*FutRet_8wk_sub2$sd[-1])
# DSpot
DSpot_8wk_sub2$sdcoef <- DSpot_8wk_sub2$fs.reg
DSpot_8wk_sub2$sdcoef$coefficients<-c(DSpot_8wk_sub2$fs.reg$coefficients[1],DSpot_8wk_sub2$fs.reg$coefficients[2:8]/DSpot_8wk_sub2$sd[1]*DSpot_8wk_sub2$sd[-1])
# DOilVol
DOilVol_8wk_sub2$sdcoef <- DOilVol_8wk_sub2$fs.reg
DOilVol_8wk_sub2$sdcoef$coefficients<-c(DOilVol_8wk_sub2$fs.reg$coefficients[1],DOilVol_8wk_sub2$fs.reg$coefficients[2:8]/DOilVol_8wk_sub2$sd[1]*DOilVol_8wk_sub2$sd[-1])
# xomRet
xomRet_8wk_sub2$sdcoef <- xomRet_8wk_sub2$fs.reg
xomRet_8wk_sub2$sdcoef$coefficients<-c(xomRet_8wk_sub2$fs.reg$coefficients[1],xomRet_8wk_sub2$fs.reg$coefficients[2:8]/xomRet_8wk_sub2$sd[1]*xomRet_8wk_sub2$sd[-1])
# bpRet
bpRet_8wk_sub2$sdcoef <- bpRet_8wk_sub2$fs.reg
bpRet_8wk_sub2$sdcoef$coefficients<-c(bpRet_8wk_sub2$fs.reg$coefficients[1],bpRet_8wk_sub2$fs.reg$coefficients[2:8]/bpRet_8wk_sub2$sd[1]*bpRet_8wk_sub2$sd[-1])
# rdsaRet
rdsaRet_8wk_sub2$sdcoef <- rdsaRet_8wk_sub2$fs.reg
rdsaRet_8wk_sub2$sdcoef$coefficients<-c(rdsaRet_8wk_sub2$fs.reg$coefficients[1],rdsaRet_8wk_sub2$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub2$sd[1]*rdsaRet_8wk_sub2$sd[-1])
# DInv
DInv_8wk_sub2$sdcoef <- DInv_8wk_sub2$fs.reg
DInv_8wk_sub2$sdcoef$coefficients<-c(DInv_8wk_sub2$fs.reg$coefficients[1],DInv_8wk_sub2$fs.reg$coefficients[2:8]/DInv_8wk_sub2$sd[1]*DInv_8wk_sub2$sd[-1])
# DProd
DProd_8wk_sub2$sdcoef <- DProd_8wk_sub2$fs.reg
DProd_8wk_sub2$sdcoef$coefficients<-c(DProd_8wk_sub2$fs.reg$coefficients[1],DProd_8wk_sub2$fs.reg$coefficients[2:8]/DProd_8wk_sub2$sd[1]*DProd_8wk_sub2$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub2$sdcoef, DSpot_8wk_sub2$sdcoef, DOilVol_8wk_sub2$sdcoef, xomRet_8wk_sub2$sdcoef, bpRet_8wk_sub2$sdcoef, rdsaRet_8wk_sub2$sdcoef, DInv_8wk_sub2$sdcoef, DProd_8wk_sub2$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub2.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub2$depVar.lagDepVar, DSpot_8wk_sub2$depVar.lagDepVar, DOilVol_8wk_sub2$depVar.lagDepVar, xomRet_8wk_sub2$depVar.lagDepVar, 
          bpRet_8wk_sub2$depVar.lagDepVar, rdsaRet_8wk_sub2$depVar.lagDepVar, DInv_8wk_sub2$depVar.lagDepVar, DProd_8wk_sub2$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub2.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub2.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub2$plot.rsq, DSpot_8wk_sub2$plot.rsq, DOilVol_8wk_sub2$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub2$plot.rsq, bpRet_8wk_sub2$plot.rsq, rdsaRet_8wk_sub2$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub2$plot.rsq, DProd_8wk_sub2$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub2.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub2$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()

#######################
#### sub3##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub3 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub3$sdcoef <- FutRet_8wk_sub3$fs.reg
FutRet_8wk_sub3$sdcoef$coefficients<-c(FutRet_8wk_sub3$fs.reg$coefficients[1],FutRet_8wk_sub3$fs.reg$coefficients[2:8]/FutRet_8wk_sub3$sd[1]*FutRet_8wk_sub3$sd[-1])
# DSpot
DSpot_8wk_sub3$sdcoef <- DSpot_8wk_sub3$fs.reg
DSpot_8wk_sub3$sdcoef$coefficients<-c(DSpot_8wk_sub3$fs.reg$coefficients[1],DSpot_8wk_sub3$fs.reg$coefficients[2:8]/DSpot_8wk_sub3$sd[1]*DSpot_8wk_sub3$sd[-1])
# DOilVol
DOilVol_8wk_sub3$sdcoef <- DOilVol_8wk_sub3$fs.reg
DOilVol_8wk_sub3$sdcoef$coefficients<-c(DOilVol_8wk_sub3$fs.reg$coefficients[1],DOilVol_8wk_sub3$fs.reg$coefficients[2:8]/DOilVol_8wk_sub3$sd[1]*DOilVol_8wk_sub3$sd[-1])
# xomRet
xomRet_8wk_sub3$sdcoef <- xomRet_8wk_sub3$fs.reg
xomRet_8wk_sub3$sdcoef$coefficients<-c(xomRet_8wk_sub3$fs.reg$coefficients[1],xomRet_8wk_sub3$fs.reg$coefficients[2:8]/xomRet_8wk_sub3$sd[1]*xomRet_8wk_sub3$sd[-1])
# bpRet
bpRet_8wk_sub3$sdcoef <- bpRet_8wk_sub3$fs.reg
bpRet_8wk_sub3$sdcoef$coefficients<-c(bpRet_8wk_sub3$fs.reg$coefficients[1],bpRet_8wk_sub3$fs.reg$coefficients[2:8]/bpRet_8wk_sub3$sd[1]*bpRet_8wk_sub3$sd[-1])
# rdsaRet
rdsaRet_8wk_sub3$sdcoef <- rdsaRet_8wk_sub3$fs.reg
rdsaRet_8wk_sub3$sdcoef$coefficients<-c(rdsaRet_8wk_sub3$fs.reg$coefficients[1],rdsaRet_8wk_sub3$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub3$sd[1]*rdsaRet_8wk_sub3$sd[-1])
# DInv
DInv_8wk_sub3$sdcoef <- DInv_8wk_sub3$fs.reg
DInv_8wk_sub3$sdcoef$coefficients<-c(DInv_8wk_sub3$fs.reg$coefficients[1],DInv_8wk_sub3$fs.reg$coefficients[2:8]/DInv_8wk_sub3$sd[1]*DInv_8wk_sub3$sd[-1])
# DProd
DProd_8wk_sub3$sdcoef <- DProd_8wk_sub3$fs.reg
DProd_8wk_sub3$sdcoef$coefficients<-c(DProd_8wk_sub3$fs.reg$coefficients[1],DProd_8wk_sub3$fs.reg$coefficients[2:8]/DProd_8wk_sub3$sd[1]*DProd_8wk_sub3$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub3$sdcoef, DSpot_8wk_sub3$sdcoef, DOilVol_8wk_sub3$sdcoef, xomRet_8wk_sub3$sdcoef, bpRet_8wk_sub3$sdcoef, rdsaRet_8wk_sub3$sdcoef, DInv_8wk_sub3$sdcoef, DProd_8wk_sub3$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub3.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub3$depVar.lagDepVar, DSpot_8wk_sub3$depVar.lagDepVar, DOilVol_8wk_sub3$depVar.lagDepVar, xomRet_8wk_sub3$depVar.lagDepVar, 
          bpRet_8wk_sub3$depVar.lagDepVar, rdsaRet_8wk_sub3$depVar.lagDepVar, DInv_8wk_sub3$depVar.lagDepVar, DProd_8wk_sub3$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub3.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub3.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub3$plot.rsq, DSpot_8wk_sub3$plot.rsq, DOilVol_8wk_sub3$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub3$plot.rsq, bpRet_8wk_sub3$plot.rsq, rdsaRet_8wk_sub3$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub3$plot.rsq, DProd_8wk_sub3$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub3.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub3$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()

#######################
#### sub4##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub4 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub4$sdcoef <- FutRet_8wk_sub4$fs.reg
FutRet_8wk_sub4$sdcoef$coefficients<-c(FutRet_8wk_sub4$fs.reg$coefficients[1],FutRet_8wk_sub4$fs.reg$coefficients[2:8]/FutRet_8wk_sub4$sd[1]*FutRet_8wk_sub4$sd[-1])
# DSpot
DSpot_8wk_sub4$sdcoef <- DSpot_8wk_sub4$fs.reg
DSpot_8wk_sub4$sdcoef$coefficients<-c(DSpot_8wk_sub4$fs.reg$coefficients[1],DSpot_8wk_sub4$fs.reg$coefficients[2:8]/DSpot_8wk_sub4$sd[1]*DSpot_8wk_sub4$sd[-1])
# DOilVol
DOilVol_8wk_sub4$sdcoef <- DOilVol_8wk_sub4$fs.reg
DOilVol_8wk_sub4$sdcoef$coefficients<-c(DOilVol_8wk_sub4$fs.reg$coefficients[1],DOilVol_8wk_sub4$fs.reg$coefficients[2:8]/DOilVol_8wk_sub4$sd[1]*DOilVol_8wk_sub4$sd[-1])
# xomRet
xomRet_8wk_sub4$sdcoef <- xomRet_8wk_sub4$fs.reg
xomRet_8wk_sub4$sdcoef$coefficients<-c(xomRet_8wk_sub4$fs.reg$coefficients[1],xomRet_8wk_sub4$fs.reg$coefficients[2:8]/xomRet_8wk_sub4$sd[1]*xomRet_8wk_sub4$sd[-1])
# bpRet
bpRet_8wk_sub4$sdcoef <- bpRet_8wk_sub4$fs.reg
bpRet_8wk_sub4$sdcoef$coefficients<-c(bpRet_8wk_sub4$fs.reg$coefficients[1],bpRet_8wk_sub4$fs.reg$coefficients[2:8]/bpRet_8wk_sub4$sd[1]*bpRet_8wk_sub4$sd[-1])
# rdsaRet
rdsaRet_8wk_sub4$sdcoef <- rdsaRet_8wk_sub4$fs.reg
rdsaRet_8wk_sub4$sdcoef$coefficients<-c(rdsaRet_8wk_sub4$fs.reg$coefficients[1],rdsaRet_8wk_sub4$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub4$sd[1]*rdsaRet_8wk_sub4$sd[-1])
# DInv
DInv_8wk_sub4$sdcoef <- DInv_8wk_sub4$fs.reg
DInv_8wk_sub4$sdcoef$coefficients<-c(DInv_8wk_sub4$fs.reg$coefficients[1],DInv_8wk_sub4$fs.reg$coefficients[2:8]/DInv_8wk_sub4$sd[1]*DInv_8wk_sub4$sd[-1])
# DProd
DProd_8wk_sub4$sdcoef <- DProd_8wk_sub4$fs.reg
DProd_8wk_sub4$sdcoef$coefficients<-c(DProd_8wk_sub4$fs.reg$coefficients[1],DProd_8wk_sub4$fs.reg$coefficients[2:8]/DProd_8wk_sub4$sd[1]*DProd_8wk_sub4$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub4$sdcoef, DSpot_8wk_sub4$sdcoef, DOilVol_8wk_sub4$sdcoef, xomRet_8wk_sub4$sdcoef, bpRet_8wk_sub4$sdcoef, rdsaRet_8wk_sub4$sdcoef, DInv_8wk_sub4$sdcoef, DProd_8wk_sub4$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub4.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub4$depVar.lagDepVar, DSpot_8wk_sub4$depVar.lagDepVar, DOilVol_8wk_sub4$depVar.lagDepVar, xomRet_8wk_sub4$depVar.lagDepVar, 
          bpRet_8wk_sub4$depVar.lagDepVar, rdsaRet_8wk_sub4$depVar.lagDepVar, DInv_8wk_sub4$depVar.lagDepVar, DProd_8wk_sub4$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub4.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub4.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub4$plot.rsq, DSpot_8wk_sub4$plot.rsq, DOilVol_8wk_sub4$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub4$plot.rsq, bpRet_8wk_sub4$plot.rsq, rdsaRet_8wk_sub4$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub4$plot.rsq, DProd_8wk_sub4$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub4.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub4$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()


#######################
#### sub5##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub5 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub5$sdcoef <- FutRet_8wk_sub5$fs.reg
FutRet_8wk_sub5$sdcoef$coefficients<-c(FutRet_8wk_sub5$fs.reg$coefficients[1],FutRet_8wk_sub5$fs.reg$coefficients[2:8]/FutRet_8wk_sub5$sd[1]*FutRet_8wk_sub5$sd[-1])
# DSpot
DSpot_8wk_sub5$sdcoef <- DSpot_8wk_sub5$fs.reg
DSpot_8wk_sub5$sdcoef$coefficients<-c(DSpot_8wk_sub5$fs.reg$coefficients[1],DSpot_8wk_sub5$fs.reg$coefficients[2:8]/DSpot_8wk_sub5$sd[1]*DSpot_8wk_sub5$sd[-1])
# DOilVol
DOilVol_8wk_sub5$sdcoef <- DOilVol_8wk_sub5$fs.reg
DOilVol_8wk_sub5$sdcoef$coefficients<-c(DOilVol_8wk_sub5$fs.reg$coefficients[1],DOilVol_8wk_sub5$fs.reg$coefficients[2:8]/DOilVol_8wk_sub5$sd[1]*DOilVol_8wk_sub5$sd[-1])
# xomRet
xomRet_8wk_sub5$sdcoef <- xomRet_8wk_sub5$fs.reg
xomRet_8wk_sub5$sdcoef$coefficients<-c(xomRet_8wk_sub5$fs.reg$coefficients[1],xomRet_8wk_sub5$fs.reg$coefficients[2:8]/xomRet_8wk_sub5$sd[1]*xomRet_8wk_sub5$sd[-1])
# bpRet
bpRet_8wk_sub5$sdcoef <- bpRet_8wk_sub5$fs.reg
bpRet_8wk_sub5$sdcoef$coefficients<-c(bpRet_8wk_sub5$fs.reg$coefficients[1],bpRet_8wk_sub5$fs.reg$coefficients[2:8]/bpRet_8wk_sub5$sd[1]*bpRet_8wk_sub5$sd[-1])
# rdsaRet
rdsaRet_8wk_sub5$sdcoef <- rdsaRet_8wk_sub5$fs.reg
rdsaRet_8wk_sub5$sdcoef$coefficients<-c(rdsaRet_8wk_sub5$fs.reg$coefficients[1],rdsaRet_8wk_sub5$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub5$sd[1]*rdsaRet_8wk_sub5$sd[-1])
# DInv
DInv_8wk_sub5$sdcoef <- DInv_8wk_sub5$fs.reg
DInv_8wk_sub5$sdcoef$coefficients<-c(DInv_8wk_sub5$fs.reg$coefficients[1],DInv_8wk_sub5$fs.reg$coefficients[2:8]/DInv_8wk_sub5$sd[1]*DInv_8wk_sub5$sd[-1])
# DProd
DProd_8wk_sub5$sdcoef <- DProd_8wk_sub5$fs.reg
DProd_8wk_sub5$sdcoef$coefficients<-c(DProd_8wk_sub5$fs.reg$coefficients[1],DProd_8wk_sub5$fs.reg$coefficients[2:8]/DProd_8wk_sub5$sd[1]*DProd_8wk_sub5$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub5$sdcoef, DSpot_8wk_sub5$sdcoef, DOilVol_8wk_sub5$sdcoef, xomRet_8wk_sub5$sdcoef, bpRet_8wk_sub5$sdcoef, rdsaRet_8wk_sub5$sdcoef, DInv_8wk_sub5$sdcoef, DProd_8wk_sub5$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub5.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub5$depVar.lagDepVar, DSpot_8wk_sub5$depVar.lagDepVar, DOilVol_8wk_sub5$depVar.lagDepVar, xomRet_8wk_sub5$depVar.lagDepVar, 
          bpRet_8wk_sub5$depVar.lagDepVar, rdsaRet_8wk_sub5$depVar.lagDepVar, DInv_8wk_sub5$depVar.lagDepVar, DProd_8wk_sub5$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub5.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub5.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub5$plot.rsq, DSpot_8wk_sub5$plot.rsq, DOilVol_8wk_sub5$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub5$plot.rsq, bpRet_8wk_sub5$plot.rsq, rdsaRet_8wk_sub5$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub5$plot.rsq, DProd_8wk_sub5$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub5.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub5$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()

#######################
#### sub6##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub6 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub6$sdcoef <- FutRet_8wk_sub6$fs.reg
FutRet_8wk_sub6$sdcoef$coefficients<-c(FutRet_8wk_sub6$fs.reg$coefficients[1],FutRet_8wk_sub6$fs.reg$coefficients[2:8]/FutRet_8wk_sub6$sd[1]*FutRet_8wk_sub6$sd[-1])
# DSpot
DSpot_8wk_sub6$sdcoef <- DSpot_8wk_sub6$fs.reg
DSpot_8wk_sub6$sdcoef$coefficients<-c(DSpot_8wk_sub6$fs.reg$coefficients[1],DSpot_8wk_sub6$fs.reg$coefficients[2:8]/DSpot_8wk_sub6$sd[1]*DSpot_8wk_sub6$sd[-1])
# DOilVol
DOilVol_8wk_sub6$sdcoef <- DOilVol_8wk_sub6$fs.reg
DOilVol_8wk_sub6$sdcoef$coefficients<-c(DOilVol_8wk_sub6$fs.reg$coefficients[1],DOilVol_8wk_sub6$fs.reg$coefficients[2:8]/DOilVol_8wk_sub6$sd[1]*DOilVol_8wk_sub6$sd[-1])
# xomRet
xomRet_8wk_sub6$sdcoef <- xomRet_8wk_sub6$fs.reg
xomRet_8wk_sub6$sdcoef$coefficients<-c(xomRet_8wk_sub6$fs.reg$coefficients[1],xomRet_8wk_sub6$fs.reg$coefficients[2:8]/xomRet_8wk_sub6$sd[1]*xomRet_8wk_sub6$sd[-1])
# bpRet
bpRet_8wk_sub6$sdcoef <- bpRet_8wk_sub6$fs.reg
bpRet_8wk_sub6$sdcoef$coefficients<-c(bpRet_8wk_sub6$fs.reg$coefficients[1],bpRet_8wk_sub6$fs.reg$coefficients[2:8]/bpRet_8wk_sub6$sd[1]*bpRet_8wk_sub6$sd[-1])
# rdsaRet
rdsaRet_8wk_sub6$sdcoef <- rdsaRet_8wk_sub6$fs.reg
rdsaRet_8wk_sub6$sdcoef$coefficients<-c(rdsaRet_8wk_sub6$fs.reg$coefficients[1],rdsaRet_8wk_sub6$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub6$sd[1]*rdsaRet_8wk_sub6$sd[-1])
# DInv
DInv_8wk_sub6$sdcoef <- DInv_8wk_sub6$fs.reg
DInv_8wk_sub6$sdcoef$coefficients<-c(DInv_8wk_sub6$fs.reg$coefficients[1],DInv_8wk_sub6$fs.reg$coefficients[2:8]/DInv_8wk_sub6$sd[1]*DInv_8wk_sub6$sd[-1])
# DProd
DProd_8wk_sub6$sdcoef <- DProd_8wk_sub6$fs.reg
DProd_8wk_sub6$sdcoef$coefficients<-c(DProd_8wk_sub6$fs.reg$coefficients[1],DProd_8wk_sub6$fs.reg$coefficients[2:8]/DProd_8wk_sub6$sd[1]*DProd_8wk_sub6$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub6$sdcoef, DSpot_8wk_sub6$sdcoef, DOilVol_8wk_sub6$sdcoef, xomRet_8wk_sub6$sdcoef, bpRet_8wk_sub6$sdcoef, rdsaRet_8wk_sub6$sdcoef, DInv_8wk_sub6$sdcoef, DProd_8wk_sub6$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub6.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub6$depVar.lagDepVar, DSpot_8wk_sub6$depVar.lagDepVar, DOilVol_8wk_sub6$depVar.lagDepVar, xomRet_8wk_sub6$depVar.lagDepVar, 
          bpRet_8wk_sub6$depVar.lagDepVar, rdsaRet_8wk_sub6$depVar.lagDepVar, DInv_8wk_sub6$depVar.lagDepVar, DProd_8wk_sub6$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub6.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub6.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub6$plot.rsq, DSpot_8wk_sub6$plot.rsq, DOilVol_8wk_sub6$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub6$plot.rsq, bpRet_8wk_sub6$plot.rsq, rdsaRet_8wk_sub6$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub6$plot.rsq, DProd_8wk_sub6$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub6.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub6$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub6$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub6$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub6$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()




#######################
#### sub7##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub7 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub7$sdcoef <- FutRet_8wk_sub7$fs.reg
FutRet_8wk_sub7$sdcoef$coefficients<-c(FutRet_8wk_sub7$fs.reg$coefficients[1],FutRet_8wk_sub7$fs.reg$coefficients[2:8]/FutRet_8wk_sub7$sd[1]*FutRet_8wk_sub7$sd[-1])
# DSpot
DSpot_8wk_sub7$sdcoef <- DSpot_8wk_sub7$fs.reg
DSpot_8wk_sub7$sdcoef$coefficients<-c(DSpot_8wk_sub7$fs.reg$coefficients[1],DSpot_8wk_sub7$fs.reg$coefficients[2:8]/DSpot_8wk_sub7$sd[1]*DSpot_8wk_sub7$sd[-1])
# DOilVol
DOilVol_8wk_sub7$sdcoef <- DOilVol_8wk_sub7$fs.reg
DOilVol_8wk_sub7$sdcoef$coefficients<-c(DOilVol_8wk_sub7$fs.reg$coefficients[1],DOilVol_8wk_sub7$fs.reg$coefficients[2:8]/DOilVol_8wk_sub7$sd[1]*DOilVol_8wk_sub7$sd[-1])
# xomRet
xomRet_8wk_sub7$sdcoef <- xomRet_8wk_sub7$fs.reg
xomRet_8wk_sub7$sdcoef$coefficients<-c(xomRet_8wk_sub7$fs.reg$coefficients[1],xomRet_8wk_sub7$fs.reg$coefficients[2:8]/xomRet_8wk_sub7$sd[1]*xomRet_8wk_sub7$sd[-1])
# bpRet
bpRet_8wk_sub7$sdcoef <- bpRet_8wk_sub7$fs.reg
bpRet_8wk_sub7$sdcoef$coefficients<-c(bpRet_8wk_sub7$fs.reg$coefficients[1],bpRet_8wk_sub7$fs.reg$coefficients[2:8]/bpRet_8wk_sub7$sd[1]*bpRet_8wk_sub7$sd[-1])
# rdsaRet
rdsaRet_8wk_sub7$sdcoef <- rdsaRet_8wk_sub7$fs.reg
rdsaRet_8wk_sub7$sdcoef$coefficients<-c(rdsaRet_8wk_sub7$fs.reg$coefficients[1],rdsaRet_8wk_sub7$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub7$sd[1]*rdsaRet_8wk_sub7$sd[-1])
# DInv
DInv_8wk_sub7$sdcoef <- DInv_8wk_sub7$fs.reg
DInv_8wk_sub7$sdcoef$coefficients<-c(DInv_8wk_sub7$fs.reg$coefficients[1],DInv_8wk_sub7$fs.reg$coefficients[2:8]/DInv_8wk_sub7$sd[1]*DInv_8wk_sub7$sd[-1])
# DProd
DProd_8wk_sub7$sdcoef <- DProd_8wk_sub7$fs.reg
DProd_8wk_sub7$sdcoef$coefficients<-c(DProd_8wk_sub7$fs.reg$coefficients[1],DProd_8wk_sub7$fs.reg$coefficients[2:8]/DProd_8wk_sub7$sd[1]*DProd_8wk_sub7$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub7$sdcoef, DSpot_8wk_sub7$sdcoef, DOilVol_8wk_sub7$sdcoef, xomRet_8wk_sub7$sdcoef, bpRet_8wk_sub7$sdcoef, rdsaRet_8wk_sub7$sdcoef, DInv_8wk_sub7$sdcoef, DProd_8wk_sub7$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub7.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub7$depVar.lagDepVar, DSpot_8wk_sub7$depVar.lagDepVar, DOilVol_8wk_sub7$depVar.lagDepVar, xomRet_8wk_sub7$depVar.lagDepVar, 
          bpRet_8wk_sub7$depVar.lagDepVar, rdsaRet_8wk_sub7$depVar.lagDepVar, DInv_8wk_sub7$depVar.lagDepVar, DProd_8wk_sub7$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub7.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub7.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub7$plot.rsq, DSpot_8wk_sub7$plot.rsq, DOilVol_8wk_sub7$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub7$plot.rsq, bpRet_8wk_sub7$plot.rsq, rdsaRet_8wk_sub7$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub7$plot.rsq, DProd_8wk_sub7$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub7.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub7$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()


#######################
#### sub8##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub8 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# --------------------------------------------------------------------------------------------------------------------------------------  
# Standardize Coefficients                                                                                  
# --------------------------------------------------------------------------------------------------------------------------------------  
##### 8wk #####
# FutRet
FutRet_8wk_sub8$sdcoef <- FutRet_8wk_sub8$fs.reg
FutRet_8wk_sub8$sdcoef$coefficients<-c(FutRet_8wk_sub8$fs.reg$coefficients[1],FutRet_8wk_sub8$fs.reg$coefficients[2:8]/FutRet_8wk_sub8$sd[1]*FutRet_8wk_sub8$sd[-1])
# DSpot
DSpot_8wk_sub8$sdcoef <- DSpot_8wk_sub8$fs.reg
DSpot_8wk_sub8$sdcoef$coefficients<-c(DSpot_8wk_sub8$fs.reg$coefficients[1],DSpot_8wk_sub8$fs.reg$coefficients[2:8]/DSpot_8wk_sub8$sd[1]*DSpot_8wk_sub8$sd[-1])
# DOilVol
DOilVol_8wk_sub8$sdcoef <- DOilVol_8wk_sub8$fs.reg
DOilVol_8wk_sub8$sdcoef$coefficients<-c(DOilVol_8wk_sub8$fs.reg$coefficients[1],DOilVol_8wk_sub8$fs.reg$coefficients[2:8]/DOilVol_8wk_sub8$sd[1]*DOilVol_8wk_sub8$sd[-1])
# xomRet
xomRet_8wk_sub8$sdcoef <- xomRet_8wk_sub8$fs.reg
xomRet_8wk_sub8$sdcoef$coefficients<-c(xomRet_8wk_sub8$fs.reg$coefficients[1],xomRet_8wk_sub8$fs.reg$coefficients[2:8]/xomRet_8wk_sub8$sd[1]*xomRet_8wk_sub8$sd[-1])
# bpRet
bpRet_8wk_sub8$sdcoef <- bpRet_8wk_sub8$fs.reg
bpRet_8wk_sub8$sdcoef$coefficients<-c(bpRet_8wk_sub8$fs.reg$coefficients[1],bpRet_8wk_sub8$fs.reg$coefficients[2:8]/bpRet_8wk_sub8$sd[1]*bpRet_8wk_sub8$sd[-1])
# rdsaRet
rdsaRet_8wk_sub8$sdcoef <- rdsaRet_8wk_sub8$fs.reg
rdsaRet_8wk_sub8$sdcoef$coefficients<-c(rdsaRet_8wk_sub8$fs.reg$coefficients[1],rdsaRet_8wk_sub8$fs.reg$coefficients[2:8]/rdsaRet_8wk_sub8$sd[1]*rdsaRet_8wk_sub8$sd[-1])
# DInv
DInv_8wk_sub8$sdcoef <- DInv_8wk_sub8$fs.reg
DInv_8wk_sub8$sdcoef$coefficients<-c(DInv_8wk_sub8$fs.reg$coefficients[1],DInv_8wk_sub8$fs.reg$coefficients[2:8]/DInv_8wk_sub8$sd[1]*DInv_8wk_sub8$sd[-1])
# DProd
DProd_8wk_sub8$sdcoef <- DProd_8wk_sub8$fs.reg
DProd_8wk_sub8$sdcoef$coefficients<-c(DProd_8wk_sub8$fs.reg$coefficients[1],DProd_8wk_sub8$fs.reg$coefficients[2:8]/DProd_8wk_sub8$sd[1]*DProd_8wk_sub8$sd[-1])

# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub8$sdcoef, DSpot_8wk_sub8$sdcoef, DOilVol_8wk_sub8$sdcoef, xomRet_8wk_sub8$sdcoef, bpRet_8wk_sub8$sdcoef, rdsaRet_8wk_sub8$sdcoef, DInv_8wk_sub8$sdcoef, DProd_8wk_sub8$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub8.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub8$depVar.lagDepVar, DSpot_8wk_sub8$depVar.lagDepVar, DOilVol_8wk_sub8$depVar.lagDepVar, xomRet_8wk_sub8$depVar.lagDepVar, 
          bpRet_8wk_sub8$depVar.lagDepVar, rdsaRet_8wk_sub8$depVar.lagDepVar, DInv_8wk_sub8$depVar.lagDepVar, DProd_8wk_sub8$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub8.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub8.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub8$plot.rsq, DSpot_8wk_sub8$plot.rsq, DOilVol_8wk_sub8$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub8$plot.rsq, bpRet_8wk_sub8$plot.rsq, rdsaRet_8wk_sub8$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub8$plot.rsq, DProd_8wk_sub8$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub8.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub8$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()


#######################
#### sub9##############
#######################
# ------------------------------------------------------------------------
# Create detrended data frame (Friday EOP)  
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
# Create detrended data frame (Tuesday EOP)     
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
           "fEpg", "sBbl", "fBbl", "sRpc", "fRpc", "sEp", "fEp", "VIX", "vix_diff", 
           "PCAsent", "PCAfreq", "PCAall", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt")


# ------------------------------------------------------------------------
# Futures Returns - 8-week horizon                                                                                         
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices, 
                        c(FutRet_t8, 
                          FutRet, DSpot, DOilVol,  OilVol, DInv, 
                          DProd, tnote_10y, DFX, sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 88302
title <- "Oil Futures Returns"

tic()
FutRet_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()



# ------------------------------------------------------------------------
# Spot Price Change - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DSpot_t8, 
                          DSpot, FutRet, DOilVol,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 933120
title <- "Spot Price Changes"

tic()
DSpot_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()

# ------------------------------------------------------------------------
# Oil Volatility - 8-week horizon
# ------------------------------------------------------------------------
subset <- dplyr::select(detrended.prices,
                        c(DOilVol_t8,
                          DOilVol, FutRet, DSpot,  OilVol, DInv, DProd, tnote_10y, DFX, 
                          sp500Ret, StkIdx, basis, WIPI_8wk, artcount, entropy, sCo, fCo, 
                          sGom, fGom, sEnv, fEnv, sEpg, fEpg, sBbl, fBbl, sRpc, 
                          fRpc, sEp, fEp, VIX, vix_diff, PCAsent, PCAfreq, PCAall,
                          BEME, Mom, BasMom, DolBeta, InflaBeta, HedgPres, liquidity, OpenInt))

horizon <- 8
sims <- 1000
seed <- 30921
title <- "Oil Volatility"

tic()
DOilVol_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 582004
title <- "Exxon Stock Returns"

tic()
xomRet_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 5285
title <- "BP Stock Returns"

tic()
bpRet_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 407783
title <- "Shell Stock Returns"

tic()
rdsaRet_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()

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

horizon <- 8
sims <- 1000
seed <- 659910
title <- "Oil Inventories"

tic()
DInv_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()


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

horizon <- 8
sims <- 1000
seed <- 7752004
title <- "Oil Production"

tic()
DProd_8wk_sub9 <- MC.FS(subset, seed, horizon, title, sims)
toc()


# --------------------------------------------------------------------------------------------------------------------------------------  
# PDF Output of All MC Simulations                                                                                    
# --------------------------------------------------------------------------------------------------------------------------------------  
# 8-week forward selection regressions with baseline, text, PCS, and risk premia measures as candidate variables
tab_model(FutRet_8wk_sub9$sdcoef, DSpot_8wk_sub9$sdcoef, DOilVol_8wk_sub9$sdcoef, xomRet_8wk_sub9$sdcoef, bpRet_8wk_sub9$sdcoef, rdsaRet_8wk_sub9$sdcoef, DInv_8wk_sub9$sdcoef, DProd_8wk_sub9$sdcoef, 
          show.ci = F, show.p = F, p.style = "star",  title = "Stepwise Forward Selection Regressions at the 8-week Horizon",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_allVars_new_sub9.doc")

# Regressions of dependent variables on lagged dependent variables 
tab_model(FutRet_8wk_sub9$depVar.lagDepVar, DSpot_8wk_sub9$depVar.lagDepVar, DOilVol_8wk_sub9$depVar.lagDepVar, xomRet_8wk_sub9$depVar.lagDepVar, 
          bpRet_8wk_sub9$depVar.lagDepVar, rdsaRet_8wk_sub9$depVar.lagDepVar, DInv_8wk_sub9$depVar.lagDepVar, DProd_8wk_sub9$depVar.lagDepVar,
          show.ci = F, show.p = F, p.style = "star",  title = "Dependent variables on lagged dependent variable coefficients",
          vcov.fun = "vcovHAC", file = "20220330/forwardSelection_8wks_AR1_new_sub9.doc")

# Adjusted R-Squared
pdf("20220330/Monte_Carlo_Exhibits_fs_rsq_new_sub9.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub9$plot.rsq, DSpot_8wk_sub9$plot.rsq, DOilVol_8wk_sub9$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub9$plot.rsq, bpRet_8wk_sub9$plot.rsq, rdsaRet_8wk_sub9$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub9$plot.rsq, DProd_8wk_sub9$plot.rsq, nrow = 3, ncol = 2)
annotate_figure(figure, top = text_grob("Adjusted R-Squared Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

dev.off()


# T-Statistics
pdf("20220330/Monte_Carlo_Exhibits_fs_tstat_new_sub9.pdf", paper = "letter", width = 7.5, height = 10, family = "Times")

figure <- ggarrange(FutRet_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit A", just = "center", face = "bold", size = 14))

figure <- ggarrange(DSpot_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit B", just = "center", face = "bold", size = 14))

figure <- ggarrange(DOilVol_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit C", just = "center", face = "bold", size = 14))

figure <- ggarrange(xomRet_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit D", just = "center", face = "bold", size = 14))

figure <- ggarrange(bpRet_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit E", just = "center", face = "bold", size = 14))

figure <- ggarrange(rdsaRet_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit F", just = "center", face = "bold", size = 14))

figure <- ggarrange(DInv_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit G", just = "center", face = "bold", size = 14))

figure <- ggarrange(DProd_8wk_sub9$plot.tstat, nrow = 2, ncol = 1)
annotate_figure(figure, top = text_grob("T-Statistic Monte Carlo Simulations:\nExhibit H", just = "center", face = "bold", size = 14))

dev.off()


#--------------------------------------------------------------------------------------- 
# Create Table A.VII "Subperiod stepwise forward selection at the eight-week horizon"
#---------------------------------------------------------------------------------------
library(rtf)
rtffile <- RTF("20220330/Table_A_VII.doc") 

#FutRet
addHeader(rtffile,title = "Table A.VII", subtitle ="Panel VII.1 : Predicting 8-weeks ahead FutRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(FutRet_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#DSpot
addHeader(rtffile,title = "", subtitle ="Panel VII.2 : Predicting 8-weeks ahead DSpot in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DSpot_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#DOilVol
addHeader(rtffile,title = "", subtitle ="Panel VII.3 : Predicting 8-weeks ahead DOilVol in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DOilVol_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#xomRet
addHeader(rtffile,title = "", subtitle ="Panel VII.4 : Predicting 8-weeks ahead xomRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(xomRet_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#bpRet
addHeader(rtffile,title = "", subtitle ="Panel VII.5 : Predicting 8-weeks ahead bpRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(bpRet_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#rdsaRet
addHeader(rtffile,title = "", subtitle ="Panel VII.6 : Predicting 8-weeks ahead rdsaRet in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(rdsaRet_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#DInv
addHeader(rtffile,title = "", subtitle ="Panel VII.7 : Predicting 8-weeks ahead DInv in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DInv_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

#DProd
addHeader(rtffile,title = "", subtitle ="Panel VII.8 : Predicting 8-weeks ahead DProd in subperiods 1-9", font.size=11)
addHeader(rtffile,title = "", subtitle ="Subperiod 1: 1998-04-01 - 1999-11-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub1$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 2: 1999-12-01 - 2002-07-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub2$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 3: 2002-08-01 - 2005-03-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub3$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 4: 2005-04-01 - 2007-11-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub4$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 5: 2007-12-01 - 2009-06-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub5$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 6: 2009-07-01 - 2012-02-29", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub6$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 7: 2012-03-01 - 2014-10-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub7$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 8: 2014-11-01 - 2017-06-30", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub8$coeff_pval), col.justify = "L")
addHeader(rtffile,title = "", subtitle ="Subperiod 9: 2017-07-01 - 2020-03-31", font.size=9)
addTable(rtffile, as.data.frame(DProd_8wk_sub9$coeff_pval), col.justify = "L")
addText(rtffile, "\n")

done(rtffile)