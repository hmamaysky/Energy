/* 
* About this .do file:
1. Transform underlying merged data as needed for regressions
2. Create 2 versions of the transformed data labeled prices or physical depending on the regressions it's meant to serve


Outputs in "./clean_data/":	
1. transformed_data_prices.dta
2. transformed_data_physical.dta
** Note that these are not the final versions of input data for the regressions because we add additional variables in subsequent .do files.
*/

********************************************************************************
* Data Transformations: prices regressions 
********************************************************************************

* ------------------- *
* Dependent variables *
* ------------------- *
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta", clear 

* Set time variable to use L. to lag vars 
tsset week 
* FutRet series are calcualted separately in update_futures_rea_v8.do.
rename FutRet_t4 FutRet_t4_Fri
rename FutRet_t8 FutRet_t8_Fri

* Oil spot price % change(we get the same result when using clc01 instead of price. They are the same series.)

gen DSpot_t4_Fri = log(price/L4.price) * 100
gen DSpot_t8_Fri = log(price/L8.price) * 100

* Exxon stock return

gen xomRet_t4_Fri = log(xomus/L4.xomus) * 100
gen xomRet_t8_Fri = log(xomus/L8.xomus) * 100

* BP stock return

gen bpRet_t4_Fri = log(bpus/L4.bpus)  * 100
gen bpRet_t8_Fri = log(bpus/L8.bpus)  * 100

* RDSA stock return 

gen rdsaRet_t4_Mon = log(rdsaus/L4.rdsaus) * 100
gen rdsaRet_t8_Mon = log(rdsaus/L8.rdsaus) * 100

* Oil price volatility difference 

gen DOilVol_t4_Fri = vol_30d - L4.vol_30d  
gen DOilVol_t8_Fri = vol_30d - L8.vol_30d  

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", replace

* --------------------- *
* Independent variables *
* --------------------- *
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", clear 

* FutRet series are calcualted separately in update_futures_rea_v8.do.
rename FutRet FutRet_Fri
rename WIPImom_4wk WIPI_4wk_monthly
rename WIPImom_8wk WIPI_8wk_monthly

* Oil price return
gen DSpot_Fri = log(price/L4.price) * 100 

* Exxon stock return
gen xomRet_Thu = log(xomus/L4.xomus) * 100

* BP stock return
gen bpRet_Thu = log(bpus/L4.bpus) * 100

* Royal Dutch Shell stock return
gen rdsaRet_Fri = log(rdsaus/L4.rdsaus) * 100

* Stock Index (Average of three energy stocks' return)
gen StkIdx_ThFr = (xomRet_Thu+bpRet_Thu+rdsaRet_Fri)/3

* S&P 500 return
gen sp500Ret_Thu = log(sp500/L4.sp500) * 100

* Oil production growth rate (We decided to follow this formula in May 2021)

gen DProd_Wed = log((prod + L.prod + L2.prod + L3.prod)/(L4.prod+ L5.prod + L6.prod+ L7.prod)) * 100

* Oil inventory growth rate
gen DInv_Wed = log(inv/L4.inv) * 100

* U.S. dollar index percent change
gen DFX_Thu = log(fx/L4.fx) * 100

* Oil price volatility first difference 
gen DOilVol_Thu = vol_30d - L4.vol_30d

* Oil price volatility level
gen OilVol_Thu = vol_30d

* VIX first difference 
*gen DVIX_Thu = vix - L4.vix

* Vix level
gen VIX_Thu = vix

* 10-year yield in levels 
* use tnote_10y
rename tnote_10y tnote_10y_Thu

* Basis 
gen basis_Fri = (clc03/clc01)^6 - 1

* Generate linear trend variables
gen trend = _n

* Text measures 4-week avg.
gen artcount_4wk = (artcount + L1.artcount + L2.artcount + L3.artcount)/4
gen entropy_4wk  = (entropy + L1.entropy + L2.entropy + L3.entropy)/4

forvalues i = 1/7 {
	gen stopic`i'_4wk  = (stopic`i' + L1.stopic`i' + L2.stopic`i' + L3.stopic`i')/4
	gen ftopic`i'_4wk  = (ftopic`i' + L1.ftopic`i' + L2.ftopic`i' + L3.ftopic`i')/4
}

* Save in excel file to use as PCA input
export excel artcount_4wk entropy_4wk stopic1_4wk stopic2_4wk stopic3_4wk stopic4_4wk stopic5_4wk stopic6_4wk stopic7_4wk using "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_prices_4wk.fri_v3.xls", firstrow(variables) replace

* Text measures: entropy
rename entropy Entropy
gen entropy_Fri = entropy_4wk 

* Text measures less entropy
* Standardize to mean 0 and unit variance 
rename artcount Artcount
egen artcount_Fri = std(artcount_4wk)
egen sCo_Fri 	  = std(stopic1_4wk)
egen fCo_Fri 	  = std(ftopic1_4wk)
egen sGom_Fri 	  = std(stopic2_4wk)
egen fGom_Fri 	  = std(ftopic2_4wk)
egen sEnv_Fri 	  = std(stopic3_4wk)
egen fEnv_Fri 	  = std(ftopic3_4wk)
egen sEpg_Fri 	  = std(stopic4_4wk)
egen fEpg_Fri 	  = std(ftopic4_4wk)
egen sBbl_Fri 	  = std(stopic5_4wk)
egen fBbl_Fri 	  = std(ftopic5_4wk)
egen sRpc_Fri 	  = std(stopic6_4wk)
egen fRpc_Fri 	  = std(ftopic6_4wk)
egen sEp_Fri 	  = std(stopic7_4wk)
egen fEp_Fri 	  = std(ftopic7_4wk)

rename (entropy_4wk artcount_4wk stopic1_4wk ftopic1_4wk stopic2_4wk ftopic2_4wk stopic3_4wk ftopic3_4wk stopic4_4wk ftopic4_4wk stopic5_4wk ftopic5_4wk stopic6_4wk ftopic6_4wk stopic7_4wk ftopic7_4wk) (entropy_4wk_Fri artcount_4wk_Fri stopic1_4wk_Fri ftopic1_4wk_Fri stopic2_4wk_Fri ftopic2_4wk_Fri stopic3_4wk_Fri ftopic3_4wk_Fri stopic4_4wk_Fri ftopic4_4wk_Fri stopic5_4wk_Fri ftopic5_4wk_Fri stopic6_4wk_Fri ftopic6_4wk_Fri stopic7_4wk_Fri ftopic7_4wk_Fri)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", nogen
drop date
rename date_wed date_Wed 
rename date_thu date_Thu 
rename date_fri date_Fri 
gen date_Tue= date_Wed-1
format date_Tue %tddd-Mon-YY
gen date_Mon= date_Fri+3
format date_Mon %tddd-Mon-YY
order date_Tue date_Wed date_Thu date_Fri date_Mon 

keep date_Tue-date_Mon monthly week FutRet_Fri-WIPI_8wk_monthly tnote_10y DSpot_Fri-fEp_Fri FutRet_t4_Fri-DOilVol_t8_Fri

order tnote_10y_Thu, after(VIX_Thu)

* Save progress so far
save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_v19.dta", replace 

* Read in PCA series created in Python to merge with our dataset.

use /Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_pca_v18, clear 

* Convert date_Fri to stata date format because date variable in python is imported as string variable in stata.
* just change date_Fri as we need just one date series to merge with transformed_data_physical_v19.dta.
gen date_Fri1 = date(date_Fri, "YMD")
format date_Fri1 %td
* drop redundant date series because we need just one date for merging.
drop date_Tue - date_Mon 
rename date_Fri1 date_Fri


save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", replace 


* Merge with transformed_data_prices_pca_v18.dta.

use /Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_v19, clear

merge 1:1 date_Fri using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", nogen

* Save progress so far
save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_v19.dta", replace 

********************************************************************************
* Data Transformations: physical regressions 
********************************************************************************

* ------------------- *
* Dependent variables *
* ------------------- *
 
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_dep.dta", clear

* Set time variable to use L. to lag vars 
tsset week 

* Oil production growth rate

gen DProd_t8_Wed  = log((prod + L.prod + L2.prod + L3.prod + L4.prod + L5.prod + L6.prod + L7.prod)/((L8.prod+ L9.prod + L10.prod+ L11.prod)*2)) * 100
gen DProd_t4_Wed  = log((prod + L.prod + L2.prod + L3.prod)/(L4.prod+ L5.prod + L6.prod+ L7.prod)) * 100


* Oil inventory growth rate 

gen DInv_t4_Wed  = log(inv/L4.inv)  * 100
gen DInv_t8_Wed  = log(inv/L8.inv)  * 100

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", replace

* --------------------- *
* Independent variables *
* --------------------- *
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", clear

* Set time variable to use L. to lag vars 
tsset week 

* FurRet is created in a separate dofile update_futures_rea_v8.do.
rename FutRet FutRet_Tue
rename WIPImom_4wk WIPI_4wk_monthly
rename WIPImom_8wk WIPI_8wk_monthly


* Oil price return
gen DSpot_Tue = log(price/L4.price) * 100 

* Exxon stock return
gen xomRet_Tue = log(xomus/L4.xomus) * 100

* BP stock return
gen bpRet_Tue = log(bpus/L4.bpus) * 100

* Royal Dutch Shell stock return
gen rdsaRet_Tue = log(rdsaus/L4.rdsaus) * 100

* Stock Index (Average of three energy stocks' return)
gen StkIdx_Tue = (xomRet_Tue+bpRet_Tue+rdsaRet_Tue)/3

* S&P 500 return
gen sp500Ret_Tue = log(sp500/L4.sp500) * 100

* Oil production growth rate (We decided to follow this formula in May 2021)

gen DProd_Wed  = log((prod + L.prod + L2.prod + L3.prod)/(L4.prod+ L5.prod + L6.prod+ L7.prod)) * 100

* Oil inventory growth rate
gen DInv_Wed = log(inv/L4.inv) * 100

* U.S. dollar index percent change
gen DFX_Tue = log(fx/L4.fx) * 100

* Oil price volatility first difference 
gen DOilVol_Tue  = vol_30d - L4.vol_30d

* Oil price volatility level
gen OilVol_Tue = vol_30d

* VIX first difference 
*gen DVIX_Tue = vix - L4.vix

* Vix level
gen VIX_Tue = vix

* 10-year yield in levels 
* use tnote_10y
rename tnote_10y tnote_10y_Tue
* Basis 
gen basis_Tue = (clc03/clc01)^6 - 1

* Generate linear trend variables
gen trend = _n

* Text measures 4-week avg.
gen artcount_4wk = (artcount + L1.artcount + L2.artcount + L3.artcount)/4
gen entropy_4wk  = (entropy + L1.entropy + L2.entropy + L3.entropy)/4

forvalues i = 1/7 {
	gen stopic`i'_4wk  = (stopic`i' + L1.stopic`i' + L2.stopic`i' + L3.stopic`i')/4
	gen ftopic`i'_4wk  = (ftopic`i' + L1.ftopic`i' + L2.ftopic`i' + L3.ftopic`i')/4
}

* Save in excel file to use as PCA input
export excel artcount_4wk entropy_4wk stopic1_4wk stopic2_4wk stopic3_4wk stopic4_4wk stopic5_4wk stopic6_4wk stopic7_4wk using "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_physical_4wk.tue_v3.xls", firstrow(variables) replace

* Text measures: entropy
rename entropy Entropy
gen entropy_Tue = entropy_4wk 

* Text measures less entropy
* Standardize to mean 0 and unit variance 
rename artcount Artcount
egen artcount_Tue = std(artcount_4wk)
egen sCo_Tue 	  = std(stopic1_4wk)
egen fCo_Tue 	  = std(ftopic1_4wk)
egen sGom_Tue 	  = std(stopic2_4wk)
egen fGom_Tue 	  = std(ftopic2_4wk)
egen sEnv_Tue 	  = std(stopic3_4wk)
egen fEnv_Tue 	  = std(ftopic3_4wk)
egen sEpg_Tue 	  = std(stopic4_4wk)
egen fEpg_Tue 	  = std(ftopic4_4wk)
egen sBbl_Tue 	  = std(stopic5_4wk)
egen fBbl_Tue 	  = std(ftopic5_4wk)
egen sRpc_Tue 	  = std(stopic6_4wk)
egen fRpc_Tue 	  = std(ftopic6_4wk)
egen sEp_Tue 	  = std(stopic7_4wk)
egen fEp_Tue 	  = std(ftopic7_4wk)

rename (entropy_4wk artcount_4wk stopic1_4wk ftopic1_4wk stopic2_4wk ftopic2_4wk stopic3_4wk ftopic3_4wk stopic4_4wk ftopic4_4wk stopic5_4wk ftopic5_4wk stopic6_4wk ftopic6_4wk stopic7_4wk ftopic7_4wk) (entropy_4wk_Tue artcount_4wk_Tue stopic1_4wk_Tue ftopic1_4wk_Tue stopic2_4wk_Tue ftopic2_4wk_Tue stopic3_4wk_Tue ftopic3_4wk_Tue stopic4_4wk_Tue ftopic4_4wk_Tue stopic5_4wk_Tue ftopic5_4wk_Tue stopic6_4wk_Tue ftopic6_4wk_Tue stopic7_4wk_Tue ftopic7_4wk_Tue)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", nogen

drop date
rename date_tue date_Tue
gen date_Wed = date_Tue+1
format date_Wed %tddd-Mon-YY

order date_Tue date_Wed


keep date_Tue date_Wed monthly week FutRet_Tue-WIPI_8wk_monthly tnote_10y_Tue DSpot_Tue-DInv_t8_Wed

drop if week>1147

* Save progress so far
save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_v19.dta", replace

* Read in PCA series created in Python to merge with our dataset.

use /Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_pca_v18, clear 

* Convert date_Tue to stata date format because date variable in python is imported as string variable in stata.
* just change date_Tue as we need just one date series to merge with transformed_data_physical_v19.dta.
gen date_Tue1 = date(date_Tue, "YMD")
format date_Tue1 %td
* drop redundant date series because we need just one date for merging.
drop date_Tue - date_Wed 
rename date_Tue1 date_Tue


save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", replace 


* Merge with transformed_data_prices_pca_v18.dta.

use /Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_v19, clear

merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/temp.dta", nogen

* Save progress so far
save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_v19.dta", replace 

