/*
This .do file creates the following weekly variables.
-Open Interest
-Liquidity 
-Hedging Pressure

Output in "./variable augmentation/":	
1. oi_tue.dta - used for RHS in physical regs
2. oi_thu.dta - used for RHS in price regs
3. liq_tue.dta - used for RHS in physical regs
4. liq_thu.dta - used for RHS in price regs
5. hp_fri.dta - used for RHS in price & physical regs

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
**************************************************************************
**For Open Interest, Hedging Pressure and Liquidity, we create it as weekly.
***************************************************************************
***** Open Interest *****
**************************
**1. Firstly, create Open Interst Thursday version.

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY


* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 4
drop dow

drop if date < daily("01/16/1992", "MDY") 
drop if date > daily("03/26/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_thu
rename price_rd price_thu
rename oiqo_rd oiqo_thu

gen week = _n

save "temp.dta", replace

***** Wednesday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 3
drop dow

drop if date < daily("01/15/1992", "MDY") 
drop if date > daily("03/25/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_wed
rename price_rd price_wed
rename oiqo_rd oiqo_wed

gen week = _n

***** Merge prices data
merge 1:1 week using "temp.dta", nogen

save "temp.dta", replace

***** Tuesday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 2
drop dow

drop if date < daily("01/14/1992", "MDY") 
drop if date > daily("03/24/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_tue
rename price_rd price_tue
rename oiqo_rd oiqo_tue

gen week = _n

***** Merge prices data
merge 1:1 week using "temp.dta", nogen


gen price = price_thu
replace price = price_wed if price == .
replace price = price_tue if price == .

gen oiqo = oiqo_thu
replace oiqo = oiqo_wed if oiqo ==.
replace oiqo = oiqo_tue if oiqo ==.

keep date_thu price oiqo

rename oiqo oiqo_thu
rename price price_thu
gen week = _n

save "temp1.dta", replace

**2. Next, we also make Open Interest Tuesday Version. 

***** Tuesday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 2
drop dow

drop if date < daily("01/14/1992", "MDY") 
drop if date > daily("03/24/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_tue
rename price_rd price_tue
rename oiqo_rd oiqo_tue

gen week = _n

save "temp2.dta", replace

***** Monday
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 1
drop dow

drop if date < daily("01/13/1992", "MDY") 
drop if date > daily("03/23/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_mon
rename price_rd price_mon
rename oiqo_rd oiqo_mon

gen week = _n

***** Merge prices data
merge 1:1 week using "temp2.dta", nogen

save "temp2.dta", replace

***** Friday(past week)
import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date
format date %tddd-Mon-YY

* Generate day of week variable to convert to weekly series
gen dow = dow(date)
keep if dow == 5
drop dow

drop if date < daily("01/10/1992", "MDY") 
drop if date > daily("03/23/2020", "MDY") 

keep date price_rd oiqo_rd
rename date date_fri
rename price_rd price_fri
rename oiqo_rd oiqo_fri

gen week = _n

***** Merge prices data
merge 1:1 week using "temp2.dta", nogen

gen price = price_tue
replace price = price_mon if price == .
replace price = price_fri if price == .

gen oiqo = oiqo_tue
replace oiqo = oiqo_mon if oiqo ==.
replace oiqo = oiqo_fri if oiqo ==.

keep date_tue price oiqo

rename oiqo oiqo_tue
rename price price_tue

gen week= _n

merge 1:1 week using "temp1.dta", nogen

* We calculate dollar open interest = OpenInt * Price.
gen oi_thu = oiqo_thu * price_thu * 1000 /*1000 is the contract size, i.e. 1 contract= 1000 barrels of oil */
gen oi_tue = oiqo_tue * price_tue * 1000

keep date_tue date_thu oi_thu oi_tue
rename date_tue date_Tue
rename date_thu date_Thu
rename oi_thu OpenInt_Thu
rename oi_tue OpenInt_Tue

save "oi_thutue.dta", replace
drop if date_Tue < daily("04/07/1998", "MDY")

drop OpenInt_Thu date_Thu
save "oi_tue",replace

use "oi_thutue.dta", clear

drop OpenInt_Tue date_Tue
save "oi_thu",replace


**********************
***** Liquiditiy *****
**********************

**1. Firstly, create Liquiditiy  Thursday version.

* ------------------------------------------------------------------------------
* FUTURES
* ------------------------------------------------------------------------------

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY

keep date clc01_rd clc02_rd vclc01_rd vclc02_rd  

* Same starting point for easy merges using week index
*drop if date < daily("04/10/1998", "MDY")

drop if date < daily("03/17/1998", "MDY")

save "clcvclc.dta", replace

* ------------------------------------------------------------------------------
* THURSDAY PRICES
* ------------------------------------------------------------------------------
gen dow = dow(date)
keep if dow == 4

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

* Thursdays signalling the expiration of the front-month contract, set dummy = 1
gen newContract_thu = 0
replace newContract_thu = 1 if day >= 21 & day <= 27
replace newContract_thu = 1 if day == 20 & month == 12
replace newContract_thu = 0 if day == 27 & month == 12

* Rename series 
rename date date_thu
rename clc01_rd clc01_thu
rename clc02_rd clc02_thu
rename vclc01_rd vclc01_thu
rename vclc02_rd vclc02_thu
rename monthly monthly_thu

order week, before(date_thu)
drop dow day month 

save "temp1.dta", replace

* ------------------------------------------------------------------------------
* WEDNESDAY PRICES
* ------------------------------------------------------------------------------
use "clcvclc.dta", clear

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
rename vclc01_rd vclc01_wed
rename vclc02_rd vclc02_wed

order week, before(date_wed)
drop dow day month monthly

merge 1:1 week using "temp1.dta", nogen

drop if week > 1150

save "temp1.dta", replace


* ------------------------------------------------------------------------------
* TUESDAY PRICES
* ------------------------------------------------------------------------------
use "clcvclc.dta", clear

gen dow = dow(date)
keep if dow == 2

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

gen newContract_tue = 0
replace newContract_tue = 1 if day >= 21 & day <= 27
replace newContract_tue = 1 if day == 20 & month == 12
replace newContract_tue = 0 if day == 27 & month == 12

* Rename series 
rename date date_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue
rename vclc01_rd vclc01_tue
rename vclc02_rd vclc02_tue

order week, before(date_tue)
drop dow day month monthly

merge 1:1 week using "temp1.dta"

drop if week > 1150

save "temp1.dta", replace

* ------------------------------------------------------------------------------
* MERGE DIFFERENT TIMINGS
* ------------------------------------------------------------------------------
* Merge series
gen clc01 = clc01_thu
replace clc01 = clc01_wed if clc01 == .
replace clc01 = clc01_tue if clc01 == .
gen clc01_dummy = 0
replace clc01_dummy = 1 if clc01_thu == . & clc01_wed != .
replace clc01_dummy = 2 if clc01_thu == . & clc01_wed == . & clc01_tue != .

gen clc02 = clc02_thu
replace clc02 = clc02_wed if clc02 == .
replace clc02 = clc02_tue if clc02 == .

gen vclc01 = vclc01_thu
replace vclc01 = vclc01_wed if vclc01 == .
replace vclc01 = vclc01_tue if vclc01 == .

gen vclc01_dummy = 0
replace vclc01_dummy = 1 if vclc01_thu == . & vclc01_wed != .
replace vclc01_dummy = 2 if vclc01_thu == . & vclc01_wed == . & vclc01_tue  != .

gen vclc02 = vclc02_thu
replace vclc02 = vclc02_wed if vclc02 == .
replace vclc02 = vclc02_tue if vclc02 == .


gen newContract = newContract_thu
replace newContract = newContract_wed if clc01_dummy == 1 
replace newContract = newContract_tue if clc01_dummy == 2 

gen newContractv = newContract_thu
replace newContractv = newContract_wed if vclc01_dummy == 1 
replace newContractv = newContract_tue if vclc01_dummy == 2

* Set time variable as week, and narrow range to 03/27/2020-3/27/2020
tsset week

drop if date_thu > daily("03/26/2020", "MDY")  
drop if date_thu < daily("03/19/1998", "MDY")   

* Compute underlying returns series
* If the front-month contract does not expire before the next Friday, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next Friday
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month


gen FutRet = (clc01 - L.clc01)/L.clc01
replace FutRet = (clc02 - L.clc02)/L.clc02 if newContract == 1 

replace FutRet = (1+ FutRet)*100

* Now we calculate the liquidity using Formula from Szymanowska et al. (2014)
replace FutRet = FutRet/100 -1

gen liquidity = vclc01/abs(FutRet)
replace liquidity = vclc02/abs(FutRet) if newContractv == 1 
/*To test whether newContarct and newContractv are the same, run the following code. 
If it is the same, we don't need to worry about any discrepancy issue.
. gen test = newContract-newContractv
. tab test
. drop test */

drop if week<4
keep date_thu liquidity
gen week=_n

rename liquidity liquidity_Thu

save "temp1.dta", replace

**2. Next, we also make Liquidity Tuesday Version. 

* ------------------------------------------------------------------------------
* FUTURES
* ------------------------------------------------------------------------------

import excel "va_raw.d.xlsx", sheet("ToStata") firstrow clear

* Convert date series from string to date

format date %tddd-Mon-YY

keep date clc01_rd clc02_rd vclc01_rd vclc02_rd  

* Same starting point for easy merges using week index
*drop if date < daily("04/10/1998", "MDY")

drop if date < daily("03/13/1998", "MDY")

save "clcvclc.dta", replace

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

* Thursdays signalling the expiration of the front-month contract, set dummy = 1
gen newContract_tue = 0
replace newContract_tue = 1 if day >= 21 & day <= 27
replace newContract_tue = 1 if day == 20 & month == 12
replace newContract_tue = 0 if day == 27 & month == 12

* Rename series 
rename date date_tue
rename clc01_rd clc01_tue
rename clc02_rd clc02_tue
rename vclc01_rd vclc01_tue
rename vclc02_rd vclc02_tue
rename monthly monthly_tue

order week, before(date_tue)
drop dow day month 

save "temp2.dta", replace

* ------------------------------------------------------------------------------
* MONDAY PRICES
* ------------------------------------------------------------------------------
use "clcvclc.dta", clear

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
rename vclc01_rd vclc01_mon
rename vclc02_rd vclc02_mon

order week, before(date_mon)
drop dow day month monthly

merge 1:1 week using "temp2.dta", nogen

drop if week > 1150

save "temp2.dta", replace

* ------------------------------------------------------------------------------
* FRIDAY PRICES (past week)
* ------------------------------------------------------------------------------
use "clcvclc.dta", clear

gen dow = dow(date)
keep if dow == 5

* Generate time variables
gen week  	= _n
gen day   	= day(date)
gen month 	= month(date)
gen monthly = mofd(date)
format monthly %tm

gen newContract_fri = 0
replace newContract_fri = 1 if day >= 21 & day <= 27
replace newContract_fri = 1 if day == 20 & month == 12
replace newContract_fri = 0 if day == 27 & month == 12

* Rename series 
rename date date_fri
rename clc01_rd clc01_fri
rename clc02_rd clc02_fri
rename vclc01_rd vclc01_fri
rename vclc02_rd vclc02_fri

order week, before(date_fri)
drop dow day month monthly

merge 1:1 week using "temp2.dta"

drop if week > 1150

save "temp2.dta", replace

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

gen vclc01 = vclc01_tue
replace vclc01 = vclc01_mon if vclc01 == .
replace vclc01 = vclc01_fri if vclc01 == .

gen vclc01_dummy = 0
replace vclc01_dummy = 1 if vclc01_tue == . & vclc01_mon != .
replace vclc01_dummy = 2 if vclc01_tue == . & vclc01_mon == . & vclc01_fri  != .

gen vclc02 = vclc02_tue
replace vclc02 = vclc02_mon if vclc02 == .
replace vclc02 = vclc02_fri if vclc02 == .


gen newContract = newContract_tue
replace newContract = newContract_mon if clc01_dummy == 1 
replace newContract = newContract_fri if clc01_dummy == 2 

gen newContractv = newContract_tue
replace newContractv = newContract_mon if vclc01_dummy == 1 
replace newContractv = newContract_fri if vclc01_dummy == 2

* Set time variable as week, and narrow range to 03/27/2020-3/27/2020
tsset week

drop if date_tue > daily("03/24/2020", "MDY")  
drop if date_tue < daily("03/17/1998", "MDY")   

* Compute underlying returns series
* If the front-month contract does not expire before the next Friday, then the return is:
* R(t) = [F1(t) - F1(t-1)]/F1(t-1), where F1 is the front-month futures series
* If the front-month contract expires before the next Friday
* R(t) = [F2(t) - F2(t-1)]/F2(t-1), where F2 is the 2-month futures series, then revert to front-month


gen FutRet = (clc01 - L.clc01)/L.clc01
replace FutRet = (clc02 - L.clc02)/L.clc02 if newContract == 1 

replace FutRet = (1+ FutRet)*100

* Now we calculate the liquidity using Formula from Szymanowska et al. (2014)
replace FutRet = FutRet/100 -1

gen liquidity = vclc01/abs(FutRet)
replace liquidity = vclc02/abs(FutRet) if newContractv == 1
drop if week<4
keep date_tue liquidity
gen week=_n

rename liquidity liquidity_Tue

merge 1:1 week using "temp1.dta", nogen


save "liq_thutue.dta", replace

drop liquidity_Thu date_thu
save "liq_tue",replace

use "liq_thutue.dta", clear

drop liquidity_Tue date_tue
save "liq_thu",replace


****************************
***** Hedging Pressure *****
****************************
*1998
import excel "cftc/Fut_1998.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT 'SWEET' - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_1998", replace

*1999
import excel "cftc/Fut_1999.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT 'SWEET' - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_1999", replace

*2000
import excel "cftc/Fut_2000.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT 'SWEET' - NEW YORK MERCANTILE EXCHANGE" | Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2000", replace

*2001
import excel "cftc/Fut_2001.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2001", replace

*2002
import excel "cftc/Fut_2002.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2002", replace

*2003
import excel "cftc/Fut_2003.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2003", replace

*2004
import excel "cftc/Fut_2004.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2004", replace

*2005
import excel "cftc/Fut_2005.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2005", replace

*2006
import excel "cftc/Fut_2006.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2006", replace

*2007
import excel "cftc/Fut_2007.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2007", replace

*2008
import excel "cftc/Fut_2008.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2008", replace

*2009
import excel "cftc/Fut_2009.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2009", replace

*2010
import excel "cftc/Fut_2010.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2010", replace

*2011
import excel "cftc/Fut_2011.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2011", replace

*2012
import excel "cftc/Fut_2012.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2012", replace

*2013
import excel "cftc/Fut_2013.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2013", replace

*2014
import excel "cftc/Fut_2014.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2014", replace

*2015
import excel "cftc/Fut_2015.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2015", replace

*2016
import excel "cftc/Fut_2016.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2016", replace

*2017
import excel "cftc/Fut_2017.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2017", replace

*2018
import excel "cftc/Fut_2018.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2018", replace

*2019
import excel "cftc/Fut_2019.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2019", replace

*2020
import excel "cftc/Fut_2020.xls", firstrow clear

keep if Market_and_Exchange_Names == "CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE"
keep Market_and_Exchange_Names Report_Date_as_MM_DD_YYYY Traders_Tot_Rept_Long_All Traders_Tot_Rept_Short_All

save "cftc/temp_2020", replace

use cftc/temp_1998.dta, clear

append using cftc/temp_1999.dta
append using cftc/temp_2000.dta
append using cftc/temp_2001.dta
append using cftc/temp_2002.dta
append using cftc/temp_2003.dta
append using cftc/temp_2004.dta
append using cftc/temp_2005.dta
append using cftc/temp_2006.dta
append using cftc/temp_2007.dta
append using cftc/temp_2008.dta
append using cftc/temp_2009.dta
append using cftc/temp_2010.dta
append using cftc/temp_2011.dta
append using cftc/temp_2012.dta
append using cftc/temp_2013.dta
append using cftc/temp_2014.dta
append using cftc/temp_2015.dta
append using cftc/temp_2016.dta
append using cftc/temp_2017.dta
append using cftc/temp_2018.dta
append using cftc/temp_2019.dta
append using cftc/temp_2020.dta

save "cftc/data_hp", replace

use cftc/data_hp, clear

drop Market_and_Exchange_Names

rename Report_Date_as_MM_DD_YYYY date
* This is Report Date, Not Data Release Date. 
* CFTC weekly data report is provided weekly on Friday at 3:30pm EST for the week ending on Tuesday
* Therefore, even though report dates are generally on Tuesdays, we regard this as a variable with Friday timing convention. 
* Because its release time is 3:30pm on Fridays, we would use the week t-1 Friday observation as the data point for both the physical and price-base week t RHS variable. 

gen dow = dow(date)

* Modify the date to match week t-1's Friday HedgPres with week t Friday/Tuesday dependent variables.

gen report_date = date /*just creating report_date to keep record of the original report date*/
format report_date %tddd-Mon-YY
 
* First, change the report dates to release dates

replace date = date+4 if dow==1
replace date = date+3 if dow==2
replace date = date+2 if dow==3
replace date = date+7 if dow ==5

* Next, do date + 7 to all observations, so that week t-1's HedgPres matches with week t's dependent variables. 

replace date = date+7

format date %tddd-Mon-YY

sort date

gen hp_fri = (Traders_Tot_Rept_Short_All-Traders_Tot_Rept_Long_All)/(Traders_Tot_Rept_Short_All+Traders_Tot_Rept_Long_All)

drop Traders_Tot_Rept_Short_All Traders_Tot_Rept_Long_All

rename hp_fri HedgPres_Fri

rename date date_Fri
gen date_Tue = date_Fri-3
format date_Tue %tddd-Mon-YY

gen week = _n
drop if week<13 | week>1158
drop week
drop dow

*It's correct that we have 1146 observation for hp. Raw data didn't have the reported observation for 9/11/2001 which was Tuesday. There was no other observation for CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE in that week. So it was not replaced. 

save hp_fri.dta, replace

