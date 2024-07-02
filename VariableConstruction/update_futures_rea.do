/* 
This .do file computes oil futures prices returns for both timing conventions. 
Code includes details on how we fill in missing observations. 
This script also includes the calculation of the World Industrial Production Index variable we use in our models.

Outputs in "./futures_rea/":												   
1. futures_rea_prices_indep.dta - used for RHS in price regs.  	   
2. futures_rea_physical_indep.dta - used for RHS in physical regs.
3. futures_rea_prices_dep - used for LHS in price regs.

*/

*-------------------------------------------------------------------------------------------------------*
* Oil Futures Prices Returns        																	* 
*-------------------------------------------------------------------------------------------------------*

***************************
**INDEPENDENT**************
* (1) Prices Regression
***************************
*We take Friday measure. 
*If Friday price is missing, we use Thurs/Wed subsequently. 

import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, clear

* Convert date series from string to date
gen date = date( v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
keep date clc01_rd clc02_rd clc03_rd

* Same starting point for easy merges using week index
drop if date < daily("04/08/1998", "MDY")

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", replace

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

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* THURSDAY PRICES
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

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

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_fri, before(date_thu)

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* WEDNESDAY PRICES
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

gen dow = dow(date)
keep if dow == 3

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
gen year	= year(date)
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

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_fri, before(date_wed)
order date_thu, before(date_wed)
drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace
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

gen clc03 = clc03_fri
replace clc03 = clc03_thu if clc03 == .
replace clc03 = clc03_wed if clc03 == .

gen newContract = newContract_fri
replace newContract = newContract_thu if clc01_dummy == 1 
replace newContract = newContract_wed if clc01_dummy == 2 

* To check that our newContract dummy = 1 at reasonable intervals comment in
* the code below:
* keep if newContract == 1
* gen time = _n
* tsset time
* gen check = week - L.week
* tab check
* Last I ran it, I got 4-5 week intervals (65.2% 4-week, and 34.8% 5-week 
* intervals between the expiration of the front-month contract)

* Set time variable as week, and narrow range to 4/10/1998-3/1/2019
tsset week
rename monthly_fri monthly

drop if date_fri > daily("03/27/2020", "MDY") 
drop if date_fri < daily("04/10/1998", "MDY")

* Compute underlying returns series
* If the front-month contract does not expire before the next Friday, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next Friday
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month
gen FutRet = (clc01 - L.clc01)/L.clc01
replace FutRet = (clc02 - L.clc02)/L.clc02 if newContract == 1 


* Independent variable: 4-week return
rename FutRet Return
gen FutRet = (1 + Return) * (1 + L.Return) * (1 + L2.Return) * (1 + L3.Return) * 100

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_indep.dta", replace

*-------------------------------------------------------------------------------------------------------*
* Baumeister and Hamilton's monthly World Industrial Production Index  		  							* 
*-------------------------------------------------------------------------------------------------------*

* Import Baumeister and Hamilton's World IP Index 
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/Baumeister_Hamilton_WIPI.xlsx", sheet("Sheet1") firstrow clear

* Generate date variables 
rename monthly date
gen monthly = mofd(date)
format monthly %tm

drop date
order monthly, before(wipi)

* Compute monthly WIPI variables of interest
* Compute the yearly change of a certain month
* Later will merge depend on lagged "monthly" variables
tsset monthly
gen WIPI    = L.wipi
gen WIPImom = (wipi - L1.wipi)/L1.wipi * 100
gen WIPIyoy = (wipi - L12.wipi)/L12.wipi * 100

* Gen monthly variable for later merging
gen monthly_4wk=monthly
format monthly_4wk %tm
gen monthly_8wk=monthly
format monthly_8wk %tm

* Drop observations outside our sample's range (note that we save the lagged months at the beginning)
drop if monthly < monthly("1998m1", "YM") | monthly > monthly("2020m3", "YM")

* Save the calculated data for later merging
save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta", replace 

*-------------------------------------------------------------------------------------------------------*
*                             Add WIPIyoy with lags into the whole dataset    	  					    * 
*-------------------------------------------------------------------------------------------------------*

* Expand monthly WIPI variables by merging with returns computed above
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_indep.dta", clear

* Generate monthly tag for 4wk analysis
gen date_4wk = date_fri-28
format date_4wk %td
gen monthly_4wks = mofd(date_4wk)
gen monthly_4wk = monthly_4wks-1
format monthly_4wk %tm
drop date_4wk
drop monthly_4wks

* Generate monthly tag for 8wk analysis
gen date_8wk = date_fri-56
format date_8wk %td
gen monthly_8wks = mofd(date_8wk)
gen monthly_8wk = monthly_8wks-1
format monthly_8wk %tm
drop date_8wk
drop monthly_8wks

* Merge WIPIyoy for 4wk analysis
merge m:1 monthly_4wk using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta"
drop wipi-WIPI
drop _merge
rename WIPImom WIPImom_4wk

* Merge WIPIyoy for 8wk analysis
merge m:1 monthly_8wk using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta"
drop wipi-WIPI
drop _merge
rename WIPImom WIPImom_8wk

sort date_fri
order date_fri, after(monthly)
order week, after(date_fri)
drop monthly_4wk-monthly_8wk

drop WIPIyoy
drop if missing(date_fri)
* Save .dta
save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_indep.dta", replace 

***************************
**INDEPENDENT**************
* (2) Phyiscal Regression
***************************
*We take Tuesday measure. 
*If Tuesday price is missing, we use prior Mon/Fri subsequently. 

import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, clear

* Convert date series from string to date
gen date = date( v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
keep date clc01_rd clc02_rd clc03_rd

* Same starting point for easy merges using week index
drop if date < daily("04/03/1998", "MDY")

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", replace

* ------------------------------------------------------------------------------
* TUESDAY PRICES
* ------------------------------------------------------------------------------
gen dow = dow(date)
keep if dow == 2

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

* TUESDAYs signalling the expiration of the front-month contract, set dummy = 1
gen newContract_tue = 0
replace newContract_tue = 1 if day >= 21 & day <= 27
replace newContract_tue = 1 if day == 20 & month == 12
replace newContract_tue = 0 if day == 27 & month == 12

* Rename series 
rename date date_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue
rename clc03_rd clc03_tue
rename monthly monthly_tue

order week, before(date_tue)
drop dow day month 

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* MONDAY PRICES
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

gen dow = dow(date)
keep if dow == 1

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
gen year	= year(date)
format monthly %tm

gen leap_year = (mod(year,4) == 0 & mod(year,100) != 0) | mod(year,400) == 0

* Mondays signalling the expiration of the front-month contract, set dummy = 1
* 	Including leap year conditions 
gen newContract_mon = 0
replace newContract_mon = 1 if day >= 24 & day <= 30

* 	Leap/common year & monday with new contract falls on 30th
replace newContract_mon = 1 if day == 2 & month == 3 & leap_year == 0
replace newContract_mon = 1 if day == 1 & month == 3 & leap_year == 1

* 	Leap/common year & Monday with new contract falls on 29th
replace newContract_mon = 1 if day == 1 & month == 3 & leap_year == 0

* 	Thanksgiving Thurs. 25 adjustment 
replace newContract_mon = 1 if day == 22 & month == 11
replace newContract_mon = 0 if day == 29 & month == 11

* 	Christmas Thurs. 25 adjustment 
replace newContract_mon = 1 if day == 22 & month == 12
replace newContract_mon = 0 if day == 29 & month == 12

* Rename series 
rename date date_mon
rename clc01_rd clc01_mon
rename clc02_rd clc02_mon
rename clc03_rd clc03_mon

order week, before(date_mon)
drop dow day month monthly year leap_year

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_tue, before(date_mon)

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* FRIDAY PRICES
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

gen dow = dow(date)
keep if dow == 5

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
gen year	= year(date)
format monthly %tm

gen newContract_fri = 0
replace newContract_fri = 1 if day >= 21 & day <= 27
replace newContract_fri = 1 if day == 20 & month == 12
replace newContract_fri = 0 if day == 27 & month == 12

* Rename series 
rename date date_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri

order week, before(date_fri)
drop dow day month monthly

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_tue, before(date_fri)
order date_mon, before(date_fri)
drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace
* ------------------------------------------------------------------------------
* MERGE DIFFERENT TIMINGS
* ------------------------------------------------------------------------------
* Merge series
gen clc01 = clc01_tue
replace clc01 = clc01_mon if clc01 == .
replace clc01 = clc01_fri if clc01 == .

gen clc01_dummy = 0
replace clc01_dummy = 1 if clc01_tue == . & clc01_mon != .
replace clc01_dummy = 2 if clc01_tue == . & clc01_mon == . & clc01_fri != .

gen clc02 = clc02_tue
replace clc02 = clc02_mon if clc02 == .
replace clc02 = clc02_fri if clc02 == .

gen clc03 = clc03_tue
replace clc03 = clc03_mon if clc02 == .
replace clc03 = clc03_fri if clc02 == .


gen newContract = newContract_tue
replace newContract = newContract_mon if clc01_dummy == 1 
replace newContract = newContract_fri if clc01_dummy == 2 

* To check that our newContract dummy = 1 at reasonable intervals comment in
* the code below:
* keep if newContract == 1
* gen time = _n
* tsset time
* gen check = week - L.week
* tab check
* Last I ran it, I got 4-5 week intervals (65.2% 4-week, and 34.8% 5-week 
* intervals between the expiration of the front-month contract)

* Set time variable as week, and narrow range to 4/10/1998-3/1/2019
tsset week
rename monthly_tue monthly

drop if date_tue > daily("03/24/2020", "MDY") 
drop if date_tue < daily("04/07/1998", "MDY")

* Compute underlying returns series
* If the front-month contract does not expire before the next TUESDAY, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next TUESDAY
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month
gen FutRet = (clc01 - L.clc01)/L.clc01
replace FutRet = (clc02 - L.clc02)/L.clc02 if newContract == 1 


* Independent variable: 4-week return
rename FutRet Return
gen FutRet = (1 + Return) * (1 + L.Return) * (1 + L2.Return) * (1 + L3.Return) * 100

drop date_mon date_fri

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_physical_indep.dta", replace

*-------------------------------------------------------------------------------------------------------*
* Baumeister and Hamilton's monthly World Industrial Production Index  		  							* 
*-------------------------------------------------------------------------------------------------------*

* Import Baumeister and Hamilton's World IP Index 
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/Baumeister_Hamilton_WIPI.xlsx", sheet("Sheet1") firstrow clear

* Generate date variables 
rename monthly date
gen monthly = mofd(date)
format monthly %tm

drop date
order monthly, before(wipi)

* Compute monthly WIPI variables of interest
* Compute the yearly change of a certain month
* Later will merge depend on lagged "monthly" variables
tsset monthly
gen WIPI    = L.wipi
gen WIPImom = (wipi - L1.wipi)/L1.wipi * 100
gen WIPIyoy = (wipi - L12.wipi)/L12.wipi * 100

* Gen monthly variable for later merging
gen monthly_4wk=monthly
format monthly_4wk %tm
gen monthly_8wk=monthly
format monthly_8wk %tm

* Drop observations outside our sample's range (note that we save the lagged months at the beginning)
drop if monthly < monthly("1998m1", "YM") | monthly > monthly("2020m3", "YM")

* Save the calculated data for later merging
save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta", replace 

*-------------------------------------------------------------------------------------------------------*
*                             Add WIPIyoy with lags into the whole dataset    	  					    * 
*-------------------------------------------------------------------------------------------------------*

* Expand monthly WIPI variables by merging with returns computed above
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_physical_indep.dta", clear

* Generate monthly tag for 4wk analysis
gen date_4wk = date_tue-28
format date_4wk %td
gen monthly_4wks = mofd(date_4wk)
gen monthly_4wk = monthly_4wks-1
format monthly_4wk %tm
drop date_4wk
drop monthly_4wks

* Generate monthly tag for 8wk analysis
gen date_8wk = date_tue-56
format date_8wk %td
gen monthly_8wks = mofd(date_8wk)
gen monthly_8wk = monthly_8wks-1
format monthly_8wk %tm
drop date_8wk
drop monthly_8wks

* Merge WIPIyoy for 4wk analysis
merge m:1 monthly_4wk using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta"
drop wipi-WIPI
drop _merge
rename WIPImom WIPImom_4wk

* Merge WIPIyoy for 8wk analysis
merge m:1 monthly_8wk using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/WIPIyoy_v2.dta"
drop wipi-WIPI
drop _merge
rename WIPImom WIPImom_8wk

sort date
order date, after(monthly)
order week, after(date)
drop monthly_4wk-monthly_8wk

drop WIPIyoy
drop if missing(date)
* Save .dta
save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_physical_indep.dta", replace


***************************
****DEPENDENT**************
*For Prices regression.
*We take Friday measure. 
*If Friday price is missing, we use next Mon/Tue subsequently. 
***************************

import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, clear

* Convert date series from string to date
gen date = date( v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
keep date clc01_rd clc02_rd

* Same starting point for easy merges using week index
drop if date < daily("04/10/1998", "MDY")

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", replace

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
rename monthly monthly_fri

order week, before(date_fri)
drop dow day month 

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* MONDAY PRICES
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

gen dow = dow(date)
keep if dow == 1

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
gen year	= year(date)
format monthly %tm

gen leap_year = (mod(year,4) == 0 & mod(year,100) != 0) | mod(year,400) == 0

* Mondays signalling the expiration of the front-month contract, set dummy = 1
* 	Including leap year conditions 
gen newContract_mon = 0
replace newContract_mon = 1 if day >= 24 & day <= 30

* 	Leap/common year & monday with new contract falls on 30th
replace newContract_mon = 1 if day == 2 & month == 3 & leap_year == 0
replace newContract_mon = 1 if day == 1 & month == 3 & leap_year == 1

* 	Leap/common year & Monday with new contract falls on 29th
replace newContract_mon = 1 if day == 1 & month == 3 & leap_year == 0

* 	Thanksgiving Thurs. 25 adjustment 
replace newContract_mon = 1 if day == 22 & month == 11
replace newContract_mon = 0 if day == 29 & month == 11

* 	Christmas Thurs. 25 adjustment 
replace newContract_mon = 1 if day == 22 & month == 12
replace newContract_mon = 0 if day == 29 & month == 12

* Rename series 
rename date date_mon
rename clc01_rd clc01_mon
rename clc02_rd clc02_mon

order week, before(date_mon)
drop dow day month monthly year leap_year

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_fri, before(date_mon)

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace

* ------------------------------------------------------------------------------
* TUESDAY PRICES - ONLY NEED THIS FOR JANUARY 4 PRICE
* THUS, NO NEED FOR NUANCED NEWCONTRACT_TUES SERIES 
* ------------------------------------------------------------------------------
use "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/daily_prices.dta", clear

gen dow = dow(date)
keep if dow == 2

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
gen year	= year(date)
format monthly %tm

* Quick newContract_tues series 
gen newContract_tue = .
replace newContract_tue = 0 if day == 4 & month == 1 & year == 2000

* Rename series 
rename date date_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue

order week, before(date_tue)
drop dow day month monthly year

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", nogen
order date_fri, before(date_tue)
order date_mon, before(date_tue)
drop if week > 1147

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/temp.dta", replace
* ------------------------------------------------------------------------------
* MERGE DIFFERENT TIMINGS
* ------------------------------------------------------------------------------
* Merge series
gen clc01 = clc01_fri
replace clc01 = clc01_mon if clc01 == .
replace clc01 = clc01_tue if clc01 == .
gen clc01_dummy = 0
replace clc01_dummy = 1 if clc01_fri == . & clc01_mon != .
replace clc01_dummy = 2 if clc01_fri == . & clc01_mon == . & clc01_tue != .

gen clc02 = clc02_fri
replace clc02 = clc02_mon if clc02 == .
replace clc02 = clc02_tue if clc02 == .

gen newContract = newContract_fri
replace newContract = newContract_mon if clc01_dummy == 1 
replace newContract = newContract_tue if clc01_dummy == 2 

* To check that our newContract dummy = 1 at reasonable intervals comment in
* the code below:
* keep if newContract == 1
* gen time = _n
* tsset time
* gen check = week - L.week
* tab check
* Last I ran it, I got 4-5 week intervals (65.2% 4-week, and 34.8% 5-week 
* intervals between the expiration of the front-month contract)

* Set time variable as week, and narrow range to 4/10/1998-3/1/2019
tsset week
rename monthly_fri monthly

drop if date_fri > daily("03/27/2020", "MDY") 
drop if date_fri < daily("04/10/1998", "MDY")

* Compute underlying returns series
* If the front-month contract does not expire before the next Friday, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next Friday
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month
gen FutRet_dep = (clc01 - L.clc01)/L.clc01
replace FutRet_dep = (clc02 - L.clc02)/L.clc02 if newContract == 1 

* Compute future returns for different horizons (4 and 8 weeks)
* 1 week ahead ( obs)

* 4 weeks ahead (829 obs)
gen FutRet_t4 = (1 + FutRet_dep) * (1 + L.FutRet_dep) * (1 + L2.FutRet_dep) * (1 + L3.FutRet_dep) * 100

* 8 weeks ahead (670 obs)
gen FutRet_t8 = (1 + FutRet_dep) * (1 + L.FutRet_dep) * (1 + L2.FutRet_dep) * (1 + L3.FutRet_dep) * (1 + L4.FutRet_dep) * ///
			    (1 + L5.FutRet_dep) * (1 + L6.FutRet_dep) * (1 + L7.FutRet_dep) * 100
				
drop date_mon date_tue

save "/Users/Economist/Dropbox/Research/ncm_research/data/futures_rea/futures_rea_prices_dep.dta", replace
		



