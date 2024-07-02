/*
This .do file creates data that is exported to python for constructing Betas(InflaBeta, DolBeta)

Output in "./variable augmentation/":	
1. cpidxyfutret_v2.dta 				 

Next, run python code "va_beta" to create InflaBeta and DolBeta

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
**Make .dta file for creating Inflation Beta and Dollar Beta. 
**We create these two in python, so we want to have .dta file with just necessary variables.
* ------------------------------------------------------------------------------

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY

****** Friday 
* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 5
drop dow

drop if date < daily("01/10/1992", "MDY") 
drop if date > daily("03/27/2020", "MDY") 


rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri
rename dxy_rd dxy_fri



keep date clc01_fri clc02_fri clc03_fri dxy_fri
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

drop if date < daily("01/09/1992", "MDY") 
drop if date > daily("03/26/2020", "MDY") 

rename date date_thu
rename clc01_rd clc01_thu
rename clc02_rd clc02_thu
rename clc03_rd clc03_thu
rename dxy_rd dxy_thu

keep date_thu clc01_thu clc02_thu clc03_thu dxy_thu
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

drop if date < daily("01/08/1992", "MDY") 
drop if date > daily("03/25/2020", "MDY") 

rename date date_wed
rename clc01_rd clc01_wed
rename clc02_rd clc02_wed
rename clc03_rd clc03_wed
rename dxy_rd dxy_wed

keep date_wed clc01_wed clc02_wed clc03_wed dxy_wed
gen week = _n

***** Merge prices data
merge 1:1 week using "temp.dta", nogen

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

keep date week clc01 clc02 clc03 dxy

gen monthly = mofd(date)
format monthly %tm


*********************************************
*Aggreate to monthly series

bysort monthly : egen  dxym = total(dxy)
bysort monthly : egen count=count(dxym)
replace dxy = dxym/count

keep monthly dxy
duplicates drop

gen dxy1 = 100*(log(dxy)-log(dxy[_n-1]))

replace dxy=dxy1

drop dxy1

merge 1:1 monthly using "futures_rea_monthly.dta", nogen
gen count=_n

save "temp.dta" , replace


import excel "cpi_yearlychg.xlsx", firstrow clear 
rename monthly m
gen count=_n

merge 1:1 count using "temp.dta", nogen
drop if missing(monthly)
drop count m FutRet2
order monthly, before(cpiyr)
save "cpidxyfutret_v2.dta", replace
