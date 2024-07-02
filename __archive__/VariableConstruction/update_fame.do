/*
This .do file reads in daily fame data saved in data/fame/fmae_raw.d.csv, and 
coverts to weekly series for both time conventions.

Outputs in "./fame/":	
1. fame_prices.thu_v4.csv - used for creating RHS in price regs.
2. fame_physical.tue_v4.csv - used for creating RHS in physical regs.
3. fame_prices.fri_v4.csv - contains the following: 
	* (1) Friday measures for constructing LHS variables in Friday timing
	* (2) Weekly prod and inv variables to construct DProd and DInv as RHS of Prices regressions

*/

*********************************************************************************************
* fame_prices.thu_v4.csv 
* RHS variable for Prices Regressions
* We use Thursday Timing Convention to avoid overlap with Friday Futures Close at 2:30pm.
* cl01 for the use of creating RHS for prices regression is calculated in update_futret_rea_v8.
*********************************************************************************************
***** Thursday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
drop v1 

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 4
drop dow price_rd clc01_rd clc02_rd clc03_rd

drop if date < daily("04/09/1998", "MDY") 
drop if date > daily("03/26/2020", "MDY") 

rename vix_rd vix_thu
rename tnote_10y_rd tnote_10y_thu
rename fx_rd fx_thu
rename sp500_rd sp500_thu

gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace

***** Wednesday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 3
drop dow price_rd clc01_rd clc02_rd clc03_rd

drop if date < daily("04/08/1998", "MDY") 
drop if date > daily("03/25/2020", "MDY") 

rename date date_wed
rename vix_rd vix_wed
rename tnote_10y_rd tnote_10y_wed
rename fx_rd fx_wed
rename sp500_rd sp500_wed

gen week = _n

***** Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace


***** Tuesday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 2
drop dow price_rd clc01_rd clc02_rd clc03_rd

drop if date < daily("04/07/1998", "MDY") 
drop if date > daily("03/24/2020", "MDY") 

rename date date_tue
rename vix_rd vix_tue
rename tnote_10y_rd tnote_10y_tue
rename fx_rd fx_tue
rename sp500_rd sp500_tue

gen week = _n

***** Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

***We fill in missing thursday values with wed/tue values subsequently.

gen vix = vix_thu
replace vix = vix_wed if vix == .
replace vix = vix_tue if vix == .

gen tnote_10y = tnote_10y_thu
replace tnote_10y = tnote_10y_wed if tnote_10y == .
replace tnote_10y = tnote_10y_tue if tnote_10y == .

gen fx = fx_thu
replace fx = fx_wed if fx == .
replace fx = fx_tue if fx == .

gen sp500 = sp500_thu
replace sp500 = sp500_wed if sp500 == .
replace sp500 = sp500_tue if sp500 == .

keep date vix tnote_10y fx sp500

export delimited using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.thu_v4.csv", replace


******************************************************************************************************************
* fame_physical.tue_v4.csv - independent variables 
* RHS for Phyiscal Regressions
* We use Tuesday Convention to avoid overlap with LHS physical variables which are released usually on Wed 10am.
******************************************************************************************************************
***** Tuesday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 2
drop dow

drop if date < daily("04/07/1998", "MDY") 
drop if date > daily("03/24/2020", "MDY") 

rename price_rd price_tue
rename vix_rd vix_tue
rename tnote_10y_rd tnote_10y_tue
rename fx_rd fx_tue
rename sp500_rd sp500_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue
rename clc03_rd clc03_tue

gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace

***** Monday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 1
drop dow

drop if date < daily("04/06/1998", "MDY") 
drop if date > daily("03/23/2020", "MDY") 

rename date date_mon
rename price_rd price_mon
rename vix_rd vix_mon
rename tnote_10y_rd tnote_10y_mon
rename fx_rd fx_mon
rename sp500_rd sp500_mon
rename clc01_rd clc01_mon
rename clc02_rd clc02_mon
rename clc03_rd clc03_mon

gen week = _n

***** Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace

***** Friday(Past week)
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 5
drop dow

drop if date < daily("04/03/1998", "MDY") 
drop if date > daily("03/20/2020", "MDY") 

rename date date_fri
rename price_rd price_fri
rename vix_rd vix_fri
rename tnote_10y_rd tnote_10y_fri
rename fx_rd fx_fri
rename sp500_rd sp500_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri

gen week = _n

***** Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

***We fill in missing tuesday values with prior mon/fri values subsequently.

gen price = price_tue
replace price = price_mon if price == .
replace price = price_fri if price == .

gen vix = vix_tue
replace vix = vix_mon if vix == .
replace vix = vix_fri if vix == .

gen tnote_10y = tnote_10y_tue
replace tnote_10y = tnote_10y_mon if tnote_10y == .
replace tnote_10y = tnote_10y_fri if tnote_10y == .

gen fx = fx_tue
replace fx = fx_mon if fx == .
replace fx = fx_fri if fx == .

gen sp500 = sp500_tue
replace sp500 = sp500_mon if sp500 == .
replace sp500 = sp500_fri if sp500 == .

gen clc01 = clc01_tue
replace clc01 = clc01_mon if clc01 == .
replace clc01 = clc01_fri if clc01 == .

gen clc02 = clc02_tue
replace clc02 = clc02_mon if clc02 == .
replace clc02 = clc02_fri if clc02 == .

gen clc03 = clc03_tue
replace clc03 = clc03_mon if clc03 == .
replace clc03 = clc03_fri if clc03 == .

keep date price vix tnote_10y fx sp500 clc01 clc02 clc03

export delimited using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_physical.tue_v4.csv", replace


***********************************************************************************************
* fame_prices.fri_v4.csv 
* this dataset contains:
* (1) Friday measures for constructing LHS variables in Friday timing
* (2) Weekly prod and inv variables to construct DProd and DInv as RHS of Prices regressions
***********************************************************************************************
***** Friday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 5
drop dow

drop if date < daily("04/10/1998", "MDY") 
drop if date > daily("03/27/2020", "MDY") 

rename price_rd price_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename clc03_rd clc03_fri

keep date price_fri clc01_fri clc02_fri clc03_fri
gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace

***** Monday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 1
drop dow

drop if date < daily("04/13/1998", "MDY") 
drop if date > daily("03/30/2020", "MDY") 

rename date date_mon
rename price_rd price_mon
rename clc01_rd clc01_mon
rename clc02_rd clc02_mon
rename clc03_rd clc03_mon

keep date price_mon clc01_mon clc02_mon clc03_mon
gen week = _n

***** Merge prices data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace

***** Tuesday
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_raw.d.csv, varnames(1) clear

* Convert date series from string to date
gen date = date(v1, "DMY", 2020)
format date %tddd-Mon-YY
order date, before(price_rd)
drop v1

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 2
drop dow

drop if date < daily("04/14/1998", "MDY") 
drop if date > daily("03/31/2020", "MDY") 

rename date date_tue
rename price_rd price_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue
rename clc03_rd clc03_tue

keep date price_tue clc01_tue clc02_tue clc03_tue
gen week = _n

***** Merge prices data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

gen price = price_fri
replace price = price_mon if price == .
replace price = price_tue if price == .

gen clc01 = clc01_fri
replace clc01 = clc01_mon if clc01 ==.
replace clc01 = clc01_tue if clc01 ==.

gen clc02 = clc02_fri
replace clc02 = clc02_mon if clc02 ==.
replace clc02 = clc02_tue if clc02 ==.

gen clc03 = clc03_fri
replace clc03 = clc03_mon if clc03 ==.
replace clc03 = clc03_tue if clc03 ==.

keep date week price clc01 clc02 clc03

save "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", replace


***** Merge prod and inv series 
import delimited /Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.fri.csv, varnames(1) clear

drop if _n > 1147

keep prod inv
gen week = _n

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/temp.dta", nogen

export delimited using "/Users/Economist/Dropbox/Research/ncm_research/data/fame/fame_prices.fri_v4.csv", replace









