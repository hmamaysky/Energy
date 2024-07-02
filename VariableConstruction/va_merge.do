/*
This .do file incorporates newly created variables in ./variable augmentation/ to the final datasets.

Output in "./clean_data/":	
1. transformed_data_physical_v19.dta
2. transformed_data_prices_v19.dta

*/

*-------------------------------------------------------------------------------------------------------*
*********************************************************************************************************
/* Setting up directory paths */

/* Sungil */
if "`c(username)'"=="Economist" {
	global folder "C:\Users\Economist\Dropbox\Research\ncm_research"
}


**********************************************************************************************************
cd "$folder"

*PHYSICAL DATA
**********************************************************************************************************
use "data\clean_data\transformed_data_physical_v19", clear 

* 1) merge new monthly variables. 
*We merge the previous month's values to the next month's weekly values.
*For example, for April 2020's weekly values, we paste March 2020's monthly value equally. 

drop monthly
gen monthly_orig = mofd(date_Tue)
format monthly_orig %tm

gen monthly= mofd(date_Tue)
format monthly %tm
replace monthly= monthly-1

save "data\clean_data\transformed_data_physical_v19", replace

merge m:1 monthly using "data\variable augmentation\data_va_monthly_v2" 
drop if _merge==2
drop monthly _merge
rename monthly_orig monthly

save "data\clean_data\transformed_data_physical_v19", replace

* 2) merge new weekly variables. 
use "data\clean_data\transformed_data_physical_v19", clear 

***hedging pressure****
merge 1:1 date_Tue using "data\variable augmentation\hp_fri", nogen

drop date_Fri report_date

save "data\clean_data\transformed_data_physical_v19", replace

***liquidity****
*merge liquidity tuesday values to phyiscal dataset. 
merge 1:1 week using "data\variable augmentation\liq_tue"

drop _merge

save "data\clean_data\transformed_data_physical_v19", replace

***open interest***
*merge open interest tuesday values to phyiscal dataset. 

merge 1:1 date_Tue using "data\variable augmentation\oi_tue"

drop _merge

**We rescale some variables
***For open interest and liquidity, we decided to take log. For momentum, we decided to subtract 100.
replace Mom_monthly = Mom_monthly-100
replace liquidity_Tue = log(liquidity_Tue)
gen OpenInt_bln_Tue= OpenInt_Tue/1000000000 /*we add this just for descriptive statistics*/
replace OpenInt_Tue = log(OpenInt_Tue)

drop date_tue

save "data\clean_data\transformed_data_physical_v19", replace


**********************************************************************************************************
*PRICE DATA
**********************************************************************************************************

use "data\clean_data\transformed_data_prices_v19", clear 

* 1) merge new monthly variables. 
*We merge the previous month's values to the next month's weekly values.
*For example, for April 2020's weekly values, we paste March 2020's monthly value equally. 

gen monthly_orig = mofd(date_Fri)
format monthly_orig %tm

replace monthly= monthly-1

merge m:1 monthly using "data\variable augmentation\data_va_monthly_v2" 
drop if _merge==2

drop monthly _merge
rename monthly_orig monthly
order monthly, before(week)

save "data\clean_data\transformed_data_prices_v19", replace


* 2) merge new weekly variables. 
use "data\clean_data\transformed_data_prices_v19", clear 

***hedging pressure(Fri)****
merge 1:1 date_Fri using "data\variable augmentation\hp_fri", nogen
drop report_date

save "data\clean_data\transformed_data_prices_v19", replace

***liquidity(Thu)****

merge 1:1 week using "data\variable augmentation\liq_thu"
drop _merge date_thu

save "data\clean_data\transformed_data_prices_v19", replace

***open interest(Thu)***
*merge open interest Thursday values to price dataset. 

merge 1:1 date_Thu using "data\variable augmentation\oi_thu"

drop _merge

***For open interest and liquidity, we decided to take log. For momentum, we decided to subtract 100.
replace Mom_monthly = Mom_monthly-100
replace liquidity_Thu = log(liquidity_Thu)
gen OpenInt_bln_Thu= OpenInt_Thu/1000000000 /*we add this just for descriptive statistics*/
replace OpenInt_Thu = log(OpenInt_Thu)

*drop redundant observations not beloning to our sample period. 
drop if missing(date_Tue)

save "data\clean_data\transformed_data_prices_v19", replace


