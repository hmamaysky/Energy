# This script counts the number of times a variable appears in a two-variable out-of-sample fixed model 
# with a run of length 3 or greater for a given dependent variable
# This script creates Table VII "Frequency of Variables Present in the Consistent Out-of-Sample Forecasting Models"

dir <- "/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/"
setwd(dir)

library(dplyr)

#FutRet

FutRet_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/FutRet.csv")

# Keep only necessary rows and columns for the table creation.
FutRet_oos_sub<- FutRet_oos_sub %>%
                  dplyr::filter(win_max_consec >=3) %>%
                  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
FutRet_oos_sub <- data.frame(Predictor=unlist(FutRet_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
FutRet_oos_sub <- FutRet_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(FutRet_oos_sub)[2] <- "FutRet"

#DSpot

DSpot_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/DSpot.csv")

# Keep only necessary rows and columns for the table creation.
DSpot_oos_sub<- DSpot_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
DSpot_oos_sub <- data.frame(Predictor=unlist(DSpot_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
DSpot_oos_sub <- DSpot_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(DSpot_oos_sub)[2] <- "DSpot"

#DOilVol

DOilVol_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/DOilVol.csv")

# Keep only necessary rows and columns for the table creation.
DOilVol_oos_sub<- DOilVol_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
DOilVol_oos_sub <- data.frame(Predictor=unlist(DOilVol_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
DOilVol_oos_sub <- DOilVol_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(DOilVol_oos_sub)[2] <- "DOilVol"

#xomRet

xomRet_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/xomRet.csv")

# Keep only necessary rows and columns for the table creation.
xomRet_oos_sub<- xomRet_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
xomRet_oos_sub <- data.frame(Predictor=unlist(xomRet_oos_sub, use.names = FALSE))

xomRet_oos_sub <- xomRet_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(xomRet_oos_sub)[2] <- "xomRet"

#bpRet

bpRet_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/bpRet.csv")

# Keep only necessary rows and columns for the table creation.
bpRet_oos_sub<- bpRet_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
bpRet_oos_sub <- data.frame(Predictor=unlist(bpRet_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
bpRet_oos_sub <- bpRet_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(bpRet_oos_sub)[2] <- "bpRet"

#rdsaRet

rdsaRet_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/rdsaRet.csv")

# Keep only necessary rows and columns for the table creation.
rdsaRet_oos_sub<- rdsaRet_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
rdsaRet_oos_sub <- data.frame(Predictor=unlist(rdsaRet_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
rdsaRet_oos_sub <- rdsaRet_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(rdsaRet_oos_sub)[2] <- "rdsaRet"

#DInv

DInv_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/DInv.csv")

# Keep only necessary rows and columns for the table creation.
DInv_oos_sub<- DInv_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
DInv_oos_sub <- data.frame(Predictor=unlist(DInv_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
DInv_oos_sub <- DInv_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(DInv_oos_sub)[2] <- "DInv"

#DProd

DProd_oos_sub <- read.csv("C:/Users/Economist/Dropbox/Research/ncm_research/OOSsubperiod/DProd.csv")

# Keep only necessary rows and columns for the table creation.
DProd_oos_sub<- DProd_oos_sub %>%
  dplyr::filter(win_max_consec >=3) %>%
  dplyr::select(var_1, var_2)

# Merge var_1 and var_2 into one column because they are irrelevant when counting.
DProd_oos_sub <- data.frame(Predictor=unlist(DProd_oos_sub, use.names = FALSE))

# Add a varaible that counts the number of appearance for each predictor
DProd_oos_sub <- DProd_oos_sub %>% add_count(Predictor) %>% distinct()

colnames(DProd_oos_sub)[2] <- "DProd"

# Merge 8 data frames for each of the dependent variable into one.

All_oos_sub <- Reduce(function(...) merge(..., all = TRUE, by = "Predictor"),
                        list(FutRet_oos_sub , DSpot_oos_sub, DOilVol_oos_sub, xomRet_oos_sub, bpRet_oos_sub, rdsaRet_oos_sub, DInv_oos_sub, DProd_oos_sub))

All_oos_sub[is.na(All_oos_sub)] <-0

#typo in variable name. StikIdx -> StkIdx
library(stringr)
All_oos_sub$Predictor =str_replace(All_oos_sub$Predictor, "StikIdx", "StkIdx")


#ord <- c(" ", "FutRet", "DSpot", "DOilVol", "xomRet", "bpRet", "rdsaRet", "StkIdx", "DInv", "DProd", "OilVol", "VIX", "DFX", "tnote_10y", "sp500Ret", 
#         "WIPI_8wk", "basis", "vix_diff", "BEME", "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity",
#         "OpenInt", "RPsdf_growing", "RPsdf_rolling", "PCAsent", "PCAfreq", "PCAall", "artcount", "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", "fEpg",
#         "sBbl", "fBbl","sRpc", "fRpc", "sEp", "fEp")
#All_oos_sub <- All_oos_sub %>% 
#  slice(match(ord, Predictor))

# We are sorting by the column 'Total' in descending order. 



All_oos_sub<- All_oos_sub %>% mutate(Total=rowSums(select_if(., is.numeric)))
All_oos_sub <-All_oos_sub[order(-All_oos_sub$Total),]


#Export the final table to .docx 

library(rtf)
rtffile <- RTF("20210827/Table_VII.doc") 
addTable(rtffile, All_oos_sub)
done(rtffile)
