# This script should run consecutively after Monte_Carlo_forwardSelection.R as it uses results created in that script.
# This script creates appendix Table A.IV with two panels 
# panel A: reports the number of selected variables of each fs regression by specification (text, nontext, all)
# panel B: reports the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)


dir <- "/Users/Economist/Dropbox/Research/ncm_research/Monte Carlo/"
setwd(dir)

# First, let's create lists of text and nontext variables for later counting purpose.

varlist_text <- c("artcount", "entropy", "sCo", "fCo", "sGom", "fGom", "sEnv", "fEnv", "sEpg", "fEpg", "sBbl", "fBbl", 
                  "sRpc", "fRpc", "sEp", "fEp", "PCAsent", "PCAfreq", "PCAall")
varlist_nontext <- c("FutRet", "DSpot", "DOilVol", "OilVol", "DInv", "xomRet", "bpRet", "rdsaRet",
                     "DProd", "tnote_10y", "DFX", "sp500Ret", "StkIdx", "basis", "WIPI_8wk","BEME", 
                     "Mom", "BasMom", "DolBeta", "InflaBeta", "HedgPres", "liquidity", "OpenInt",
                     "VIX", "vix_diff")

#FutRet

pval <-list(FutRet_8wk$plot.tstat$plot_env$pval1,FutRet_8wk$plot.tstat$plot_env$pval2,
            FutRet_8wk$plot.tstat$plot_env$pval3,FutRet_8wk$plot.tstat$plot_env$pval4,
            FutRet_8wk$plot.tstat$plot_env$pval5,FutRet_8wk$plot.tstat$plot_env$pval6,
            FutRet_8wk$plot.tstat$plot_env$pval7,FutRet_8wk$plot.tstat$plot_env$pval8)

pval_FutRet <- data.table(Predictor=names(FutRet_8wk$fs.reg$coefficients), Pval=pval)

pval_FutRet <- pval_FutRet[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_FutRet <- pval_FutRet %>%
  mutate(text = ifelse(pval_FutRet$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_FutRet$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_FutRet$Predictor %in% varlist_text & pval_FutRet$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_FutRet$Predictor %in% varlist_nontext & pval_FutRet$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#DSpot

pval <-list(DSpot_8wk$plot.tstat$plot_env$pval1,DSpot_8wk$plot.tstat$plot_env$pval2,
            DSpot_8wk$plot.tstat$plot_env$pval3,DSpot_8wk$plot.tstat$plot_env$pval4,
            DSpot_8wk$plot.tstat$plot_env$pval5,DSpot_8wk$plot.tstat$plot_env$pval6,
            DSpot_8wk$plot.tstat$plot_env$pval7,DSpot_8wk$plot.tstat$plot_env$pval8)

pval_DSpot <- data.table(Predictor=names(DSpot_8wk$fs.reg$coefficients), Pval=pval)

pval_DSpot <- pval_DSpot[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_DSpot <- pval_DSpot %>%
  mutate(text = ifelse(pval_DSpot$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_DSpot$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_DSpot$Predictor %in% varlist_text & pval_DSpot$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_DSpot$Predictor %in% varlist_nontext & pval_DSpot$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#DOilVol

pval <-list(DOilVol_8wk$plot.tstat$plot_env$pval1,DOilVol_8wk$plot.tstat$plot_env$pval2,
            DOilVol_8wk$plot.tstat$plot_env$pval3,DOilVol_8wk$plot.tstat$plot_env$pval4,
            DOilVol_8wk$plot.tstat$plot_env$pval5,DOilVol_8wk$plot.tstat$plot_env$pval6,
            DOilVol_8wk$plot.tstat$plot_env$pval7,DOilVol_8wk$plot.tstat$plot_env$pval8)

pval_DOilVol <- data.table(Predictor=names(DOilVol_8wk$fs.reg$coefficients), Pval=pval)

pval_DOilVol <- pval_DOilVol[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_DOilVol <- pval_DOilVol %>%
  mutate(text = ifelse(pval_DOilVol$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_DOilVol$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_DOilVol$Predictor %in% varlist_text & pval_DOilVol$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_DOilVol$Predictor %in% varlist_nontext & pval_DOilVol$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#xomRet

pval <-list(xomRet_8wk$plot.tstat$plot_env$pval1,xomRet_8wk$plot.tstat$plot_env$pval2,
            xomRet_8wk$plot.tstat$plot_env$pval3,xomRet_8wk$plot.tstat$plot_env$pval4,
            xomRet_8wk$plot.tstat$plot_env$pval5,xomRet_8wk$plot.tstat$plot_env$pval6,
            xomRet_8wk$plot.tstat$plot_env$pval7,xomRet_8wk$plot.tstat$plot_env$pval8)

pval_xomRet <- data.table(Predictor=names(xomRet_8wk$fs.reg$coefficients), Pval=pval)

pval_xomRet <- pval_xomRet[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_xomRet <- pval_xomRet %>%
  mutate(text = ifelse(pval_xomRet$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_xomRet$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_xomRet$Predictor %in% varlist_text & pval_xomRet$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_xomRet$Predictor %in% varlist_nontext & pval_xomRet$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#bpRet

pval <-list(bpRet_8wk$plot.tstat$plot_env$pval1,bpRet_8wk$plot.tstat$plot_env$pval2,
            bpRet_8wk$plot.tstat$plot_env$pval3,bpRet_8wk$plot.tstat$plot_env$pval4,
            bpRet_8wk$plot.tstat$plot_env$pval5,bpRet_8wk$plot.tstat$plot_env$pval6,
            bpRet_8wk$plot.tstat$plot_env$pval7,bpRet_8wk$plot.tstat$plot_env$pval8)

pval_bpRet <- data.table(Predictor=names(bpRet_8wk$fs.reg$coefficients), Pval=pval)

pval_bpRet <- pval_bpRet[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_bpRet <- pval_bpRet %>%
  mutate(text = ifelse(pval_bpRet$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_bpRet$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_bpRet$Predictor %in% varlist_text & pval_bpRet$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_bpRet$Predictor %in% varlist_nontext & pval_bpRet$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#rdsaRet

pval <-list(rdsaRet_8wk$plot.tstat$plot_env$pval1,rdsaRet_8wk$plot.tstat$plot_env$pval2,
            rdsaRet_8wk$plot.tstat$plot_env$pval3,rdsaRet_8wk$plot.tstat$plot_env$pval4,
            rdsaRet_8wk$plot.tstat$plot_env$pval5,rdsaRet_8wk$plot.tstat$plot_env$pval6,
            rdsaRet_8wk$plot.tstat$plot_env$pval7,rdsaRet_8wk$plot.tstat$plot_env$pval8)

pval_rdsaRet <- data.table(Predictor=names(rdsaRet_8wk$fs.reg$coefficients), Pval=pval)

pval_rdsaRet <- pval_rdsaRet[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_rdsaRet <- pval_rdsaRet %>%
  mutate(text = ifelse(pval_rdsaRet$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_rdsaRet$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_rdsaRet$Predictor %in% varlist_text & pval_rdsaRet$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_rdsaRet$Predictor %in% varlist_nontext & pval_rdsaRet$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#DInv

pval <-list(DInv_8wk$plot.tstat$plot_env$pval1,DInv_8wk$plot.tstat$plot_env$pval2,
            DInv_8wk$plot.tstat$plot_env$pval3,DInv_8wk$plot.tstat$plot_env$pval4,
            DInv_8wk$plot.tstat$plot_env$pval5,DInv_8wk$plot.tstat$plot_env$pval6,
            DInv_8wk$plot.tstat$plot_env$pval7,DInv_8wk$plot.tstat$plot_env$pval8)

pval_DInv <- data.table(Predictor=names(DInv_8wk$fs.reg$coefficients), Pval=pval)

pval_DInv <- pval_DInv[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_DInv <- pval_DInv %>%
  mutate(text = ifelse(pval_DInv$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_DInv$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_DInv$Predictor %in% varlist_text & pval_DInv$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_DInv$Predictor %in% varlist_nontext & pval_DInv$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)

#DProd

pval <-list(DProd_8wk$plot.tstat$plot_env$pval1,DProd_8wk$plot.tstat$plot_env$pval2,
            DProd_8wk$plot.tstat$plot_env$pval3,DProd_8wk$plot.tstat$plot_env$pval4,
            DProd_8wk$plot.tstat$plot_env$pval5,DProd_8wk$plot.tstat$plot_env$pval6,
            DProd_8wk$plot.tstat$plot_env$pval7,DProd_8wk$plot.tstat$plot_env$pval8)

pval_DProd <- data.table(Predictor=names(DProd_8wk$fs.reg$coefficients), Pval=pval)

pval_DProd <- pval_DProd[Predictor != '(Intercept)']

# Count the following: 
# 1. the number of selected variables of each fs regression by specification (text, nontext, all)
# 2. the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_DProd <- pval_DProd %>%
  mutate(text = ifelse(pval_DProd$Predictor %in% varlist_text, 1, 0),
         nontext= ifelse(pval_DProd$Predictor %in% varlist_nontext, 1, 0),
         text_signif = ifelse(pval_DProd$Predictor %in% varlist_text & pval_DProd$Pval <0.05, 1, 0),
         nontext_signif = ifelse(pval_DProd$Predictor %in% varlist_nontext & pval_DProd$Pval <0.05, 1, 0)
  ) %>%
  summarize(count_text = sum(text),
            count_nontext = sum(nontext),
            count_all = count_text+count_nontext,
            count_text_sig = sum(text_signif),
            count_nontext_sig = sum(nontext_signif),
            count_all_sig = count_text_sig+count_nontext_sig)


# now create the table
# panel A: reports the number of selected variables of each fs regression by specification (text, nontext, all)


count_selected <- data.table(Specification=c("text","nontext","all"), FutRet=as.numeric(as.data.frame(count_FutRet)[1,1:3]),
                             DSpot=as.numeric(as.data.frame(count_DSpot)[1,1:3]),
                             DOilVol=as.numeric(as.data.frame(count_DOilVol)[1,1:3]),
                             xomRet=as.numeric(as.data.frame(count_xomRet)[1,1:3]),
                             bpRet=as.numeric(as.data.frame(count_bpRet)[1,1:3]),
                             rdsaRet=as.numeric(as.data.frame(count_rdsaRet)[1,1:3]),
                             DInv=as.numeric(as.data.frame(count_DInv)[1,1:3]),
                             DProd=as.numeric(as.data.frame(count_DProd)[1,1:3]))


# Add a column 'All' for rowsum.
count_selected <- left_join(count_selected, count_selected[,2:9] %>%
                              mutate(All = rowSums(.)))


# panel B: reports the number of selected variables of each fs regression that are significant, by specification (text, nontext, all)
count_selected_sig <- data.table(Specification=c("text","nontext","all"), FutRet=as.numeric(as.data.frame(count_FutRet)[1,4:6]),
                             DSpot=as.numeric(as.data.frame(count_DSpot)[1,4:6]),
                             DOilVol=as.numeric(as.data.frame(count_DOilVol)[1,4:6]),
                             xomRet=as.numeric(as.data.frame(count_xomRet)[1,4:6]),
                             bpRet=as.numeric(as.data.frame(count_bpRet)[1,4:6]),
                             rdsaRet=as.numeric(as.data.frame(count_rdsaRet)[1,4:6]),
                             DInv=as.numeric(as.data.frame(count_DInv)[1,4:6]),
                             DProd=as.numeric(as.data.frame(count_DProd)[1,4:6]))

# Add a column 'All' for rowsum.
count_selected_sig <- left_join(count_selected_sig, count_selected_sig[,2:9] %>%
                                  mutate(All = rowSums(.)))


#Export the final table to .docx 

library(rtf)
rtffile <- RTF("20210821/Table_A_IV.doc") 

addHeader(rtffile,title = "Table A. IV ", subtitle ="Count of Selected and Statistically Significant Variables of Forward Selection Model", font.size=11)
addHeader(rtffile,title = "", subtitle ="Panel A: Count of selected variables", font.size=9)
addTable(rtffile, as.data.frame(count_selected))
addHeader(rtffile,title = "", subtitle ="Panel B: Count of selected and significant Variables", font.size=9)
addTable(rtffile, as.data.frame(count_selected_sig))

done(rtffile)

