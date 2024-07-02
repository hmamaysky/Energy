/*
About this .do file:
1. I create a "week" variable to identify weekly observations 
since within the same week obs. fall on diff. days for diff. vars
2. I create four merged datasets, two following the time conventions
used in oil price returns, oil company stock returns, and oil volatility regressions, 
and two others for oil prod. and inv.
3. The datasets are in ./clean_data/ labeled for prices and physical
4. The data in the merged_data .dta files is the underlying data for all transformations, 
thus, only transformations need to change when considering different time horizons

Outputs in "./clean_data/":	
1. merged_data_prices_dep.dta
2. merged_data_prices_indep.dta
3. merged_data_physical_dep.dta
4. merged_data_physical_indep.data
*/

********************************************************************************
* Merging all relevant data: 
* (1) dataset for LHS  in Prices regressions 
********************************************************************************

* import and format oil prices data
* -----------------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.fri_v4.csv, varnames(1) clear 

*gen week = _n 
*drop v1
drop date prod inv

order week, before(price)

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta", replace 


* import and format stock price and oil price volatility data
* -----------------------------------------------------------
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_dep.fri.xls", sheet("Sheet1") firstrow clear

gen week = _n 
format date %td 

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta"

order week, before(xomus)
drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta", replace 


* import and format royal dutch shell stock price data
* -----------------------------------------------------------

import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_dep.mon.xls", sheet("Sheet1") firstrow clear

gen week = _n 
format date %td 
drop date 

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta"

order week, before(xomus)
order rdsaus, after(bpus)
drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta", replace 

* import and format futures returns
* -----------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_dep.dta", clear

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta"

drop _merge

drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_dep.dta", replace


********************************************************************************
* Merging all relevant data: 
* (2) dataset for RHS  in prices regressions 
********************************************************************************

* import and format oil prices, inv and prod data
* -----------------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.fri_v4.csv, varnames(1) clear 

*gen week = _n 
*drop v1
drop date
keep week prod inv

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", replace 


* import and format stock price and oil price volatility data
* -----------------------------------------------------------
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_indep.frithu.xls", sheet("Sheet1") firstrow clear

gen week = _n 
format date %td 

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta"

order week, before(xomus)
drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", replace 


* import and format financial variables 
* -------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.thu_v4.csv, varnames(1) clear 

gen week = _n 
*drop v1
drop date
order week, before(vix)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta"

drop _merge
order date, before(week)

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", replace 


* import and format text measures 
* -------------------------------
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_prices.fri_v3.xlsx", sheet("Sheet1") firstrow clear

gen week = _n 
drop date
order week, before(artcount)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta"

drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", replace 


* import and format futures returns and WIPI data
* -----------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_indep.dta", clear

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta"

drop _merge

drop if week > 1147

gen price = clc01 /*not really necessary because we use clc01 for both Futures Price and Spot Price, but just create a seaparte series for convenience and consistency. They are exactly the same.*/

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_prices_indep.dta", replace


********************************************************************************
* Merging all relevant data: 
* (3) dataset for LHS  in physical regressions 
********************************************************************************

* import and format oil inv and prod data
* ---------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_physical.fri.csv, varnames(1) clear 

gen week = _n 
gen date = date( v1, "DMY", 2020)
format date %tddd-Mon-YY
drop v1

order week, before(prod)

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_dep.dta", replace 

********************************************************************************
* Merging all relevant data: 
* (4) dataset for RHS  in physical regressions 
********************************************************************************

* import and format oil inv and prod data
* ---------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_physical.fri.csv, varnames(1) clear 

gen week = _n 
gen date = date( v1, "DMY", 2020)
format date %tddd-Mon-YY
drop v1

order week, before(prod)

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace 

* import and format stock price data 
* -----------------------------------
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_physical_indep.tue.xls", sheet("Sheet1") firstrow clear

gen week = _n
drop date

merge 1:1 week using  "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta"

drop if week > 1147
drop _merge 

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace 

* import and format financial variables 
* -------------------------------------
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_physical.tue_v4.csv, varnames(1) clear 

gen week = _n 
*drop v1
drop date
order week, before(price)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta"

drop _merge
order date, before(week)

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace 

* import and format text measures 
* -------------------------------
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_physical.tue_v3.xlsx", sheet("Sheet1") firstrow clear

gen week = _n 
drop date
order week, before(artcount)

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta"

drop _merge
drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace 

* import and format futures and WIPI data
* --------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_physical_indep.dta", clear

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta"

drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace

drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/merged_data_physical_indep.dta", replace