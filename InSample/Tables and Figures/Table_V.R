# This script is to measure the instability of In-Sample forward selection model. 
# This script creates Table V "Measuring Instability of In-Sample Forward Selection Results"
# Should run consecutively after stepwiseForwardSelection_detrended_sub because this script uses the results from it.

library(dplyr)
library(data.table)

### FutRet

# Read in data
FutRet_sub1 <- data.table(Predictor=names(FutRet_t8_sub1$stdcoef),Coeff_sub1=(FutRet_t8_sub1$stdcoef))
FutRet_sub2 <- data.table(Predictor=names(FutRet_t8_sub2$stdcoef),Coeff_sub2=(FutRet_t8_sub2$stdcoef))
FutRet_sub3 <- data.table(Predictor=names(FutRet_t8_sub3$stdcoef),Coeff_sub3=(FutRet_t8_sub3$stdcoef))
FutRet_sub4 <- data.table(Predictor=names(FutRet_t8_sub4$stdcoef),Coeff_sub4=(FutRet_t8_sub4$stdcoef))
FutRet_sub5 <- data.table(Predictor=names(FutRet_t8_sub5$stdcoef),Coeff_sub5=(FutRet_t8_sub5$stdcoef))
FutRet_sub6 <- data.table(Predictor=names(FutRet_t8_sub6$stdcoef),Coeff_sub6=(FutRet_t8_sub6$stdcoef))
FutRet_sub7 <- data.table(Predictor=names(FutRet_t8_sub7$stdcoef),Coeff_sub7=(FutRet_t8_sub7$stdcoef))
FutRet_sub8 <- data.table(Predictor=names(FutRet_t8_sub8$stdcoef),Coeff_sub8=(FutRet_t8_sub8$stdcoef))
FutRet_sub9 <- data.table(Predictor=names(FutRet_t8_sub9$stdcoef),Coeff_sub9=(FutRet_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
FutRet_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
       list(FutRet_sub1, FutRet_sub2 , FutRet_sub3, FutRet_sub4, FutRet_sub5, FutRet_sub6, FutRet_sub7, FutRet_sub8, FutRet_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
FutRet_sub_all <- merge(FutRet_sub_all,var_list_full, all = TRUE, by = "Predictor")
FutRet_sub_all <- FutRet_sub_all[Predictor != '(Intercept)']
FutRet_sub_all[is.na(FutRet_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_FutRet <- cor(FutRet_sub_all[,!c("Predictor")])
lowtri_FutRet <- cormat_FutRet[lower.tri(cormat_FutRet)]
avgcor_FutRet <- round(mean(lowtri_FutRet),2)

# Count positive and negative coefficients separately.
count_pos_FutRet <- melt(FutRet_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_FutRet, "V1", "pos")
count_neg_FutRet <- melt(FutRet_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_FutRet, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_FutRet <- data.table(Predictor=FutRet_sub_all$Predictor, FutRet=paste0(count_pos_FutRet$pos, ",", count_neg_FutRet$neg))

# Add average correlation as the last row of the table
count_posneg_corr_FutRet <- rbind(count_posneg_FutRet,list("Avg. corr.", avgcor_FutRet))


### DSpot

# Read in data
DSpot_sub1 <- data.table(Predictor=names(DSpot_t8_sub1$stdcoef),Coeff_sub1=(DSpot_t8_sub1$stdcoef))
DSpot_sub2 <- data.table(Predictor=names(DSpot_t8_sub2$stdcoef),Coeff_sub2=(DSpot_t8_sub2$stdcoef))
DSpot_sub3 <- data.table(Predictor=names(DSpot_t8_sub3$stdcoef),Coeff_sub3=(DSpot_t8_sub3$stdcoef))
DSpot_sub4 <- data.table(Predictor=names(DSpot_t8_sub4$stdcoef),Coeff_sub4=(DSpot_t8_sub4$stdcoef))
DSpot_sub5 <- data.table(Predictor=names(DSpot_t8_sub5$stdcoef),Coeff_sub5=(DSpot_t8_sub5$stdcoef))
DSpot_sub6 <- data.table(Predictor=names(DSpot_t8_sub6$stdcoef),Coeff_sub6=(DSpot_t8_sub6$stdcoef))
DSpot_sub7 <- data.table(Predictor=names(DSpot_t8_sub7$stdcoef),Coeff_sub7=(DSpot_t8_sub7$stdcoef))
DSpot_sub8 <- data.table(Predictor=names(DSpot_t8_sub8$stdcoef),Coeff_sub8=(DSpot_t8_sub8$stdcoef))
DSpot_sub9 <- data.table(Predictor=names(DSpot_t8_sub9$stdcoef),Coeff_sub9=(DSpot_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
DSpot_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                        list(DSpot_sub1, DSpot_sub2 , DSpot_sub3, DSpot_sub4, DSpot_sub5, DSpot_sub6, DSpot_sub7, DSpot_sub8, DSpot_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
DSpot_sub_all <- merge(DSpot_sub_all,var_list_full, all = TRUE, by = "Predictor")
DSpot_sub_all = DSpot_sub_all[Predictor != '(Intercept)']
DSpot_sub_all[is.na(DSpot_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_DSpot <- cor(DSpot_sub_all[,!c("Predictor")])
lowtri_DSpot <- cormat_DSpot[lower.tri(cormat_DSpot)]
avgcor_DSpot <- round(mean(lowtri_DSpot),2)

# Count positive and negative coefficients separately.
count_pos_DSpot <- melt(DSpot_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_DSpot, "V1", "pos")
count_neg_DSpot <- melt(DSpot_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_DSpot, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_DSpot <- data.table(Predictor=DSpot_sub_all$Predictor, DSpot=paste0(count_pos_DSpot$pos, ",", count_neg_DSpot$neg))

# Add average correlation as the last row of the table
count_posneg_corr_DSpot <- rbind(count_posneg_DSpot,list("Avg. corr.", avgcor_DSpot))

### DOilVol

# Read in data
DOilVol_sub1 <- data.table(Predictor=names(DOilVol_t8_sub1$stdcoef),Coeff_sub1=(DOilVol_t8_sub1$stdcoef))
DOilVol_sub2 <- data.table(Predictor=names(DOilVol_t8_sub2$stdcoef),Coeff_sub2=(DOilVol_t8_sub2$stdcoef))
DOilVol_sub3 <- data.table(Predictor=names(DOilVol_t8_sub3$stdcoef),Coeff_sub3=(DOilVol_t8_sub3$stdcoef))
DOilVol_sub4 <- data.table(Predictor=names(DOilVol_t8_sub4$stdcoef),Coeff_sub4=(DOilVol_t8_sub4$stdcoef))
DOilVol_sub5 <- data.table(Predictor=names(DOilVol_t8_sub5$stdcoef),Coeff_sub5=(DOilVol_t8_sub5$stdcoef))
DOilVol_sub6 <- data.table(Predictor=names(DOilVol_t8_sub6$stdcoef),Coeff_sub6=(DOilVol_t8_sub6$stdcoef))
DOilVol_sub7 <- data.table(Predictor=names(DOilVol_t8_sub7$stdcoef),Coeff_sub7=(DOilVol_t8_sub7$stdcoef))
DOilVol_sub8 <- data.table(Predictor=names(DOilVol_t8_sub8$stdcoef),Coeff_sub8=(DOilVol_t8_sub8$stdcoef))
DOilVol_sub9 <- data.table(Predictor=names(DOilVol_t8_sub9$stdcoef),Coeff_sub9=(DOilVol_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
DOilVol_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                          list(DOilVol_sub1, DOilVol_sub2 , DOilVol_sub3, DOilVol_sub4, DOilVol_sub5, DOilVol_sub6, DOilVol_sub7, DOilVol_sub8, DOilVol_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
DOilVol_sub_all <- merge(DOilVol_sub_all,var_list_full, all = TRUE, by = "Predictor")
DOilVol_sub_all = DOilVol_sub_all[Predictor != '(Intercept)']
DOilVol_sub_all[is.na(DOilVol_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_DOilVol <- cor(DOilVol_sub_all[,!c("Predictor")])
lowtri_DOilVol <- cormat_DOilVol[lower.tri(cormat_DOilVol)]
avgcor_DOilVol <- round(mean(lowtri_DOilVol),2)

# Count positive and negative coefficients separately.
count_pos_DOilVol <- melt(DOilVol_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_DOilVol, "V1", "pos")
count_neg_DOilVol <- melt(DOilVol_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_DOilVol, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_DOilVol <- data.table(Predictor=DOilVol_sub_all$Predictor, DOilVol=paste0(count_pos_DOilVol$pos, ",", count_neg_DOilVol$neg))

# Add average correlation as the last row of the table
count_posneg_corr_DOilVol <- rbind(count_posneg_DOilVol,list("Avg. corr.", avgcor_DOilVol))


### xomRet

# Read in data
xomRet_sub1 <- data.table(Predictor=names(xomRet_t8_sub1$stdcoef),Coeff_sub1=(xomRet_t8_sub1$stdcoef))
xomRet_sub2 <- data.table(Predictor=names(xomRet_t8_sub2$stdcoef),Coeff_sub2=(xomRet_t8_sub2$stdcoef))
xomRet_sub3 <- data.table(Predictor=names(xomRet_t8_sub3$stdcoef),Coeff_sub3=(xomRet_t8_sub3$stdcoef))
xomRet_sub4 <- data.table(Predictor=names(xomRet_t8_sub4$stdcoef),Coeff_sub4=(xomRet_t8_sub4$stdcoef))
xomRet_sub5 <- data.table(Predictor=names(xomRet_t8_sub5$stdcoef),Coeff_sub5=(xomRet_t8_sub5$stdcoef))
xomRet_sub6 <- data.table(Predictor=names(xomRet_t8_sub6$stdcoef),Coeff_sub6=(xomRet_t8_sub6$stdcoef))
xomRet_sub7 <- data.table(Predictor=names(xomRet_t8_sub7$stdcoef),Coeff_sub7=(xomRet_t8_sub7$stdcoef))
xomRet_sub8 <- data.table(Predictor=names(xomRet_t8_sub8$stdcoef),Coeff_sub8=(xomRet_t8_sub8$stdcoef))
xomRet_sub9 <- data.table(Predictor=names(xomRet_t8_sub9$stdcoef),Coeff_sub9=(xomRet_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
xomRet_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                         list(xomRet_sub1, xomRet_sub2 , xomRet_sub3, xomRet_sub4, xomRet_sub5, xomRet_sub6, xomRet_sub7, xomRet_sub8, xomRet_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
xomRet_sub_all <- merge(xomRet_sub_all,var_list_full, all = TRUE, by = "Predictor")
xomRet_sub_all = xomRet_sub_all[Predictor != '(Intercept)']
xomRet_sub_all[is.na(xomRet_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_xomRet <- cor(xomRet_sub_all[,!c("Predictor")])
lowtri_xomRet <- cormat_xomRet[lower.tri(cormat_xomRet)]
avgcor_xomRet <- round(mean(lowtri_xomRet),2)

# Count positive and negative coefficients separately.
count_pos_xomRet <- melt(xomRet_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_xomRet, "V1", "pos")
count_neg_xomRet <- melt(xomRet_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_xomRet, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_xomRet <- data.table(Predictor=xomRet_sub_all$Predictor, xomRet=paste0(count_pos_xomRet$pos, ",", count_neg_xomRet$neg))

# Add average correlation as the last row of the table
count_posneg_corr_xomRet <- rbind(count_posneg_xomRet,list("Avg. corr.", avgcor_xomRet))



### bpRet

# Read in data
bpRet_sub1 <- data.table(Predictor=names(bpRet_t8_sub1$stdcoef),Coeff_sub1=(bpRet_t8_sub1$stdcoef))
bpRet_sub2 <- data.table(Predictor=names(bpRet_t8_sub2$stdcoef),Coeff_sub2=(bpRet_t8_sub2$stdcoef))
bpRet_sub3 <- data.table(Predictor=names(bpRet_t8_sub3$stdcoef),Coeff_sub3=(bpRet_t8_sub3$stdcoef))
bpRet_sub4 <- data.table(Predictor=names(bpRet_t8_sub4$stdcoef),Coeff_sub4=(bpRet_t8_sub4$stdcoef))
bpRet_sub5 <- data.table(Predictor=names(bpRet_t8_sub5$stdcoef),Coeff_sub5=(bpRet_t8_sub5$stdcoef))
bpRet_sub6 <- data.table(Predictor=names(bpRet_t8_sub6$stdcoef),Coeff_sub6=(bpRet_t8_sub6$stdcoef))
bpRet_sub7 <- data.table(Predictor=names(bpRet_t8_sub7$stdcoef),Coeff_sub7=(bpRet_t8_sub7$stdcoef))
bpRet_sub8 <- data.table(Predictor=names(bpRet_t8_sub8$stdcoef),Coeff_sub8=(bpRet_t8_sub8$stdcoef))
bpRet_sub9 <- data.table(Predictor=names(bpRet_t8_sub9$stdcoef),Coeff_sub9=(bpRet_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
bpRet_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                        list(bpRet_sub1, bpRet_sub2 , bpRet_sub3, bpRet_sub4, bpRet_sub5, bpRet_sub6, bpRet_sub7, bpRet_sub8, bpRet_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
bpRet_sub_all <- merge(bpRet_sub_all,var_list_full, all = TRUE, by = "Predictor")
bpRet_sub_all = bpRet_sub_all[Predictor != '(Intercept)']
bpRet_sub_all[is.na(bpRet_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_bpRet <- cor(bpRet_sub_all[,!c("Predictor")])
lowtri_bpRet <- cormat_bpRet[lower.tri(cormat_bpRet)]
avgcor_bpRet <- round(mean(lowtri_bpRet),2)

# Count positive and negative coefficients separately.
count_pos_bpRet <- melt(bpRet_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_bpRet, "V1", "pos")
count_neg_bpRet <- melt(bpRet_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_bpRet, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_bpRet <- data.table(Predictor=bpRet_sub_all$Predictor, bpRet=paste0(count_pos_bpRet$pos, ",", count_neg_bpRet$neg))

# Add average correlation as the last row of the table
count_posneg_corr_bpRet <- rbind(count_posneg_bpRet,list("Avg. corr.", avgcor_bpRet))


### rdsaRet

# Read in data
rdsaRet_sub1 <- data.table(Predictor=names(rdsaRet_t8_sub1$stdcoef),Coeff_sub1=(rdsaRet_t8_sub1$stdcoef))
rdsaRet_sub2 <- data.table(Predictor=names(rdsaRet_t8_sub2$stdcoef),Coeff_sub2=(rdsaRet_t8_sub2$stdcoef))
rdsaRet_sub3 <- data.table(Predictor=names(rdsaRet_t8_sub3$stdcoef),Coeff_sub3=(rdsaRet_t8_sub3$stdcoef))
rdsaRet_sub4 <- data.table(Predictor=names(rdsaRet_t8_sub4$stdcoef),Coeff_sub4=(rdsaRet_t8_sub4$stdcoef))
rdsaRet_sub5 <- data.table(Predictor=names(rdsaRet_t8_sub5$stdcoef),Coeff_sub5=(rdsaRet_t8_sub5$stdcoef))
rdsaRet_sub6 <- data.table(Predictor=names(rdsaRet_t8_sub6$stdcoef),Coeff_sub6=(rdsaRet_t8_sub6$stdcoef))
rdsaRet_sub7 <- data.table(Predictor=names(rdsaRet_t8_sub7$stdcoef),Coeff_sub7=(rdsaRet_t8_sub7$stdcoef))
rdsaRet_sub8 <- data.table(Predictor=names(rdsaRet_t8_sub8$stdcoef),Coeff_sub8=(rdsaRet_t8_sub8$stdcoef))
rdsaRet_sub9 <- data.table(Predictor=names(rdsaRet_t8_sub9$stdcoef),Coeff_sub9=(rdsaRet_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
rdsaRet_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                          list(rdsaRet_sub1, rdsaRet_sub2 , rdsaRet_sub3, rdsaRet_sub4, rdsaRet_sub5, rdsaRet_sub6, rdsaRet_sub7, rdsaRet_sub8, rdsaRet_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.prices)[c(13:17, 21:26, 28, 30:58)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
rdsaRet_sub_all <- merge(rdsaRet_sub_all,var_list_full, all = TRUE, by = "Predictor")
rdsaRet_sub_all = rdsaRet_sub_all[Predictor != '(Intercept)']
rdsaRet_sub_all[is.na(rdsaRet_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_rdsaRet <- cor(rdsaRet_sub_all[,!c("Predictor")])
lowtri_rdsaRet <- cormat_rdsaRet[lower.tri(cormat_rdsaRet)]
avgcor_rdsaRet <- round(mean(lowtri_rdsaRet),2)

# Count positive and negative coefficients separately.
count_pos_rdsaRet <- melt(rdsaRet_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_rdsaRet, "V1", "pos")
count_neg_rdsaRet <- melt(rdsaRet_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_rdsaRet, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_rdsaRet <- data.table(Predictor=rdsaRet_sub_all$Predictor, rdsaRet=paste0(count_pos_rdsaRet$pos, ",", count_neg_rdsaRet$neg))

# Add average correlation as the last row of the table
count_posneg_corr_rdsaRet <- rbind(count_posneg_rdsaRet,list("Avg. corr.", avgcor_rdsaRet))



### DInv

# Read in data
DInv_sub1 <- data.table(Predictor=names(DInv_t8_sub1$stdcoef),Coeff_sub1=(DInv_t8_sub1$stdcoef))
DInv_sub2 <- data.table(Predictor=names(DInv_t8_sub2$stdcoef),Coeff_sub2=(DInv_t8_sub2$stdcoef))
DInv_sub3 <- data.table(Predictor=names(DInv_t8_sub3$stdcoef),Coeff_sub3=(DInv_t8_sub3$stdcoef))
DInv_sub4 <- data.table(Predictor=names(DInv_t8_sub4$stdcoef),Coeff_sub4=(DInv_t8_sub4$stdcoef))
DInv_sub5 <- data.table(Predictor=names(DInv_t8_sub5$stdcoef),Coeff_sub5=(DInv_t8_sub5$stdcoef))
DInv_sub6 <- data.table(Predictor=names(DInv_t8_sub6$stdcoef),Coeff_sub6=(DInv_t8_sub6$stdcoef))
DInv_sub7 <- data.table(Predictor=names(DInv_t8_sub7$stdcoef),Coeff_sub7=(DInv_t8_sub7$stdcoef))
DInv_sub8 <- data.table(Predictor=names(DInv_t8_sub8$stdcoef),Coeff_sub8=(DInv_t8_sub8$stdcoef))
DInv_sub9 <- data.table(Predictor=names(DInv_t8_sub9$stdcoef),Coeff_sub9=(DInv_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
DInv_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                       list(DInv_sub1, DInv_sub2 , DInv_sub3, DInv_sub4, DInv_sub5, DInv_sub6, DInv_sub7, DInv_sub8, DInv_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.physical)[c(5:7, 11:18, 20, 22:50)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
DInv_sub_all <- merge(DInv_sub_all,var_list_full, all = TRUE, by = "Predictor")
DInv_sub_all = DInv_sub_all[Predictor != '(Intercept)']
DInv_sub_all[is.na(DInv_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_DInv <- cor(DInv_sub_all[,!c("Predictor")])
lowtri_DInv <- cormat_DInv[lower.tri(cormat_DInv)]
avgcor_DInv <- round(mean(lowtri_DInv),2)

# Count positive and negative coefficients separately.
count_pos_DInv <- melt(DInv_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_DInv, "V1", "pos")
count_neg_DInv <- melt(DInv_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_DInv, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_DInv <- data.table(Predictor=DInv_sub_all$Predictor, DInv=paste0(count_pos_DInv$pos, ",", count_neg_DInv$neg))

# Add average correlation as the last row of the table
count_posneg_corr_DInv <- rbind(count_posneg_DInv,list("Avg. corr.", avgcor_DInv))


### DProd

# Read in data
DProd_sub1 <- data.table(Predictor=names(DProd_t8_sub1$stdcoef),Coeff_sub1=(DProd_t8_sub1$stdcoef))
DProd_sub2 <- data.table(Predictor=names(DProd_t8_sub2$stdcoef),Coeff_sub2=(DProd_t8_sub2$stdcoef))
DProd_sub3 <- data.table(Predictor=names(DProd_t8_sub3$stdcoef),Coeff_sub3=(DProd_t8_sub3$stdcoef))
DProd_sub4 <- data.table(Predictor=names(DProd_t8_sub4$stdcoef),Coeff_sub4=(DProd_t8_sub4$stdcoef))
DProd_sub5 <- data.table(Predictor=names(DProd_t8_sub5$stdcoef),Coeff_sub5=(DProd_t8_sub5$stdcoef))
DProd_sub6 <- data.table(Predictor=names(DProd_t8_sub6$stdcoef),Coeff_sub6=(DProd_t8_sub6$stdcoef))
DProd_sub7 <- data.table(Predictor=names(DProd_t8_sub7$stdcoef),Coeff_sub7=(DProd_t8_sub7$stdcoef))
DProd_sub8 <- data.table(Predictor=names(DProd_t8_sub8$stdcoef),Coeff_sub8=(DProd_t8_sub8$stdcoef))
DProd_sub9 <- data.table(Predictor=names(DProd_t8_sub9$stdcoef),Coeff_sub9=(DProd_t8_sub9$stdcoef))

# Merge data from 9 subperiods into a single data table
DProd_sub_all <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                        list(DProd_sub1, DProd_sub2 , DProd_sub3, DProd_sub4, DProd_sub5, DProd_sub6, DProd_sub7, DProd_sub8, DProd_sub9))

# Elaborate list of variables that are relevant for the subperiod forward selection 
var_list_full <- colnames(subset.physical)[c(5:7, 11:18, 20, 22:50)]
var_list_full <- data.table(Predictor = var_list_full)

# This step is conducted to show the full list of variables even though it was never seleected in any of 9 subperiods
DProd_sub_all <- merge(DProd_sub_all,var_list_full, all = TRUE, by = "Predictor")
DProd_sub_all = DProd_sub_all[Predictor != '(Intercept)']
DProd_sub_all[is.na(DProd_sub_all)] <-0

# Calculate average correlation of the 9 subperiod vectors of coefficients. 
cormat_DProd <- cor(DProd_sub_all[,!c("Predictor")])
lowtri_DProd <- cormat_DProd[lower.tri(cormat_DProd)]
avgcor_DProd <- round(mean(lowtri_DProd),2)

# Count positive and negative coefficients separately.
count_pos_DProd <- melt(DProd_sub_all, id.vars="Predictor")[,sum(value>0), by=Predictor]
setnames(count_pos_DProd, "V1", "pos")
count_neg_DProd <- melt(DProd_sub_all, id.vars="Predictor")[,sum(value<0), by=Predictor]
setnames(count_neg_DProd, "V1", "neg")

# Make a table for ordered pairs of the count of pos/neg coefficients
count_posneg_DProd <- data.table(Predictor=DProd_sub_all$Predictor, DProd=paste0(count_pos_DProd$pos, ",", count_neg_DProd$neg))

# Add average correlation as the last row of the table
count_posneg_corr_DProd <- rbind(count_posneg_DProd,list("Avg. corr.", avgcor_DProd))


# Merge 8 data tables created above for each LHS.
# We create a single table that includes all the count information and avg. corr. information for the 8 lhs varaibles.

List_DT=list(count_posneg_corr_FutRet,count_posneg_corr_DSpot,count_posneg_corr_DOilVol,count_posneg_corr_xomRet,count_posneg_corr_bpRet,count_posneg_corr_rdsaRet,count_posneg_corr_DInv,count_posneg_corr_DProd)

MergedDT <- Reduce(function(...) merge(..., all = TRUE), List_DT)

MergedDT <- rbind(MergedDT,list(" ", "+,-", "+,-", "+,-", "+,-", "+,-", "+,-", "+,-", "+,-"))

ord <- c(" ", "FutRet", "DSpot", "DOilVol", "StkIdx", "DInv", "DProd", "OilVol", "VIX", "DFX", "tnote_10y", "sp500Ret", 
           "WIPI_8wk", "basis", "vix_diff", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity",
           "OpenInt", "PCAsent", "PCAfreq", "PCAall", "artcount", "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", "fEpg",
           "sBbl", "fBbl","sRpc", "fRpc", "sEp", "fEp", "Avg. corr.")
MergedDT <- MergedDT %>% 
            slice(match(ord, Predictor))

MergedDT[MergedDT == "0,0"] <- " "

#Export the final table to .docx 

library(rtf)
rtffile <- RTF("/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/20210821/Table_V.doc") 
addTable(rtffile, as.data.frame(MergedDT))
done(rtffile)

