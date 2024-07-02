*-------------------------------------------------------------------------------------------------------*
* Variable Augmentation - monthly variables				 
/*
This .do file creates the following monthly variables.
-BE/ME
-Momentum
-Basis Momentum

Output in "./variable augmentation/":	
1. data_va_monthly_v2.dta

*/

*-------------------------------------------------------------------------------------------------------*
*********************************************************************************************************
/* Setting up directory paths */

/* Sungil */
if "`c(username)'"=="Economist" {
	global folder "C:\Users\Economist\Dropbox\Research\ncm_research"
}


**********************************************************************************************************
cd "$folder\data\variable augmentation"

* ------------------------------------------------------------------------------
* FUTURES
* ------------------------------------------------------------------------------

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY

keep date clc01_rd clc02_rd clc03_rd

* Same starting point for easy merges using week index
*drop if date < daily("04/10/1998", "MDY")

drop if date < daily("01/08/1992", "MDY")

save "daily_prices_new.dta", replace

* ------------------------------------------------------------------------------
* FRIDAY PRICES
* ------------------------------------------------------------------------------
gen dow = dow(date)
keep if dow == 5

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

* Fridays signalling the expiration of the front-month contract, set dummy = 1
gen newContract_fri = 0
replace newContract_fri = 1 if day >= 21 & day <= 27
replace newContract_fri = 1 if day == 20 & month == 12
replace newContract_fri = 0 if day == 27 & month == 12

* Rename series 
rename date date_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri
rename monthly monthly_fri

order week, before(date_fri)
drop dow day month 

save "temp_new.dta", replace

* ------------------------------------------------------------------------------
* THURSDAY PRICES
* ------------------------------------------------------------------------------
use "daily_prices_new.dta", clear

gen dow = dow(date)
keep if dow == 4

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

gen newContract_thu = 0
replace newContract_thu = 1 if day >= 21 & day <= 27
replace newContract_thu = 1 if day == 20 & month == 12
replace newContract_thu = 0 if day == 27 & month == 12

* Rename series 
rename date date_thu
rename clc01_rd clc01_thu
rename clc02_rd clc02_thu
rename clc03_rd clc03_thu

order week, before(date_thu)
drop dow day month monthly

merge 1:1 week using "temp_new.dta", nogen
order date_fri, before(date_thu)

save "temp_new.dta", replace

* ------------------------------------------------------------------------------
* WEDNESDAY PRICES
* ------------------------------------------------------------------------------

use "daily_prices_new.dta", clear

gen dow = dow(date)
keep if dow == 3

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

gen newContract_wed = 0
replace newContract_wed = 1 if day >= 21 & day <= 27
replace newContract_wed = 1 if day == 20 & month == 12
replace newContract_wed = 0 if day == 27 & month == 12

* Rename series 
rename date date_wed
rename clc01_rd clc01_wed
rename clc02_rd clc02_wed
rename clc03_rd clc03_wed

order week, before(date_wed)
drop dow day month monthly

merge 1:1 week using "temp_new.dta", nogen
order date_fri, before(date_wed)
order date_thu, before(date_wed)
*drop if week > 1147
save "temp_new.dta", replace
* ------------------------------------------------------------------------------
* MERGE DIFFERENT TIMINGS
* ------------------------------------------------------------------------------
* Merge series
gen clc01 = clc01_fri
replace clc01 = clc01_thu if clc01 == .
replace clc01 = clc01_wed if clc01 == .
gen clc01_dummy = 0
replace clc01_dummy = 1 if clc01_fri == . & clc01_thu != .
replace clc01_dummy = 2 if clc01_fri == . & clc01_thu == . & clc01_wed != .

gen clc02 = clc02_fri
replace clc02 = clc02_thu if clc02 == .
replace clc02 = clc02_wed if clc02 == .

*gen clc02_dummy = 0
*replace clc02_dummy = 1 if clc02_fri == . & clc02_thu != .
*replace clc02_dummy = 2 if clc02_fri == . & clc02_thu == . & clc02_wed != .
*do not need to create clc02_dummy since it's just the same as clc01_dummy

gen clc03 = clc03_fri
replace clc03 = clc03_thu if clc03 == .
replace clc03 = clc03_wed if clc03 == .

gen newContract = newContract_fri
replace newContract = newContract_thu if clc01_dummy == 1 
replace newContract = newContract_wed if clc01_dummy == 2 


* Set time variable as week, and narrow range
tsset week
rename date_fri date
rename monthly_fri monthly

drop if date > daily("03/27/2020", "MDY")  
drop if date < daily("01/08/1992", "MDY")   

* Compute underlying returns series
* If the front-month contract does not expire before the next Friday, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next Friday
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month


gen FutRet1 = (clc01 - L.clc01)/L.clc01
replace FutRet1 = (clc02 - L.clc02)/L.clc02 if newContract == 1 

gen FutRet2 = (clc02 - L.clc02)/L.clc02
replace FutRet2 = (clc03 - L.clc03)/L.clc03 if newContract == 1 

gen FutRet1_1 = 1+ FutRet1
gen FutRet2_1 = 1+ FutRet2


drop if missing(FutRet1_1) | missing(FutRet2_1)
replace FutRet1_1=log(FutRet1_1)
replace FutRet2_1=log(FutRet2_1)


keep monthly week date clc01 clc02 clc03 FutRet1_1 FutRet2_1

*weekly to monthly : we do this for calculating new variables like BM, Momentum, Basis Momentum which need monthly variables.
bysort monthly : egen  product1 = total(FutRet1_1)
replace product1 = exp(product1)
bysort monthly : egen  product2 = total(FutRet2_1)
replace product2 = exp(product2)

gen FutRet1 = (product1) * 100
gen FutRet2 = (product2) * 100
 
keep monthly FutRet1 FutRet2
duplicates drop
save "futures_rea_monthly.dta", replace 

* ------------------------------------------------------------------------------
* SPOT AND OTHER QUANTITY VARIABLES
* ------------------------------------------------------------------------------

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY
***** Friday
* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 5
drop dow

drop if date > daily("03/27/2020", "MDY")  
drop if date < daily("01/10/1992", "MDY")  

rename price_rd price_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri
rename dxy_rd dxy_fri
rename oiqo_rd oiqo_fri


keep date price_fri clc01_fri clc02_fri clc03_fri dxy_fri oiqo_fri
gen week = _n

save "temp.dta", replace

***** Thursday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 4
drop dow

drop if date > daily("03/26/2020", "MDY")  
drop if date < daily("01/09/1992", "MDY")  

rename date date_thu
rename price_rd price_thu
rename clc01_rd clc01_thu
rename clc02_rd clc02_thu
rename clc03_rd clc03_thu
rename dxy_rd dxy_thu
rename oiqo_rd oiqo_thu

keep date price_thu clc01_thu clc02_thu clc03_thu dxy_thu oiqo_thu
gen week = _n

***** Merge prices data
merge 1:1 week using "temp.dta", nogen

save "temp.dta", replace

***** Wednesday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 3
drop dow

drop if date > daily("03/25/2020", "MDY")  
drop if date < daily("01/08/1992", "MDY")  

rename date date_wed
rename price_rd price_wed
rename clc01_rd clc01_wed
rename clc02_rd clc02_wed
rename clc03_rd clc03_wed
rename dxy_rd dxy_wed
rename oiqo_rd oiqo_wed

keep date price_wed clc01_wed clc02_wed clc03_wed dxy_wed oiqo_wed
gen week = _n

***** Merge prices data
merge 1:1 week using "temp.dta", nogen

gen price = price_fri
replace price = price_thu if price == .
replace price = price_wed if price == .

gen clc01 = clc01_fri
replace clc01 = clc01_thu if clc01 ==.
replace clc01 = clc01_wed if clc01 ==.

gen clc02 = clc02_fri
replace clc02 = clc02_thu if clc02 ==.
replace clc02 = clc02_wed if clc02 ==.

gen clc03 = clc03_fri
replace clc03 = clc03_thu if clc03 ==.
replace clc03 = clc03_wed if clc03 ==.

gen dxy = dxy_fri
replace dxy = dxy_thu if dxy ==.
replace dxy = dxy_wed if dxy ==.

gen oiqo = oiqo_fri
replace oiqo = oiqo_thu if oiqo ==.
replace oiqo = oiqo_wed if oiqo ==.

keep date week price clc01 clc02 clc03 dxy oiqo
gen monthly = mofd(date)
format monthly %tm

save "temp.dta", replace


**Aggregate price from weekly to monthly( averaging)
bysort monthly : egen sum2= total(price)
bysort monthly : egen count=count(sum2)
replace price= sum2/count

keep monthly price
duplicates drop

merge 1:1 monthly using "futures_rea_monthly.dta", nogen
order monthly, before(price)


* ------------------------------------------------------------------------------
* NEW VARIABLE CONSTRUCTION
* ------------------------------------------------------------------------------
***************
***** BE/ME*******
***************
gen count=_n
xtset price count
order count, before(monthly)
rangestat (mean) wanted=price, interval(count -67 -56)

sort count

rename wanted avg12
gen bmratio = avg12/price[_n-1] /*Stata cannot name a variable with "/" character. 
So we need to name it BE/ME manually later for Tables and Figures*/

drop avg12

***************
***Momentum****
***************

*We create momentum just for Futures Returns.
*Following Boons et al. (2017) We use the common measure of the prior 12-month cumulative raw return on the asset, skipping the most recent month’s return. (from month t-11 to t-1)
*****Futures
replace FutRet1 = FutRet1/100
gen mom_fut1 = 100* FutRet1[_n-11] * FutRet1[_n-10] * FutRet1[_n-9] * FutRet1[_n-8] * FutRet1[_n-7] * FutRet1[_n-6] * FutRet1[_n-5] * FutRet1[_n-4] * FutRet1[_n-3] * FutRet1[_n-2] * FutRet1[_n-1]
replace FutRet1 = 100* FutRet1

****************
*Basis Momentum*
****************

replace FutRet2 = FutRet2/100
 
gen mom_fut2 = 100*FutRet2[_n-11] * FutRet2[_n-10] * FutRet2[_n-9] * FutRet2[_n-8] * FutRet2[_n-7] * FutRet2[_n-6] * FutRet2[_n-5] * FutRet2[_n-4] * FutRet2[_n-3] * FutRet2[_n-2] * FutRet2[_n-1]
replace FutRet2 = 100* FutRet2

gen basismom= mom_fut1 - mom_fut2

rename bmratio BEME_monthly /*The variable name we show in the tables is BE/ME, but Stata doesn't allow "/" in variable names. We adjust the name later manually*/
rename mom_fut1 Mom_monthly
rename basismom BasMom_monthly

keep count monthly BEME_monthly Mom_monthly BasMom_monthly

save "data_va_monthly_v2.dta", replace


***Merge Beta variables created in python

merge 1:1 count using "dxy_beta_v2.dta", nogen

drop index

save "data_va_monthly_v2.dta", replace

merge 1:1 count using "cpi_beta_v2.dta", nogen

drop index

*Keep only from March 1998
drop if count<75
drop count 

rename dxy_betas DolBeta_monthly
rename cpiyr_betas InflaBeta_monthly

save "data_va_monthly_v2.dta", replace
