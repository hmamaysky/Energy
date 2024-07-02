/* 
This .do file reads in daily Bloomberg data saved in data/bloomberg/bloomberg raw.d.xlsx, 
and creates follwoing four output files with weekly data for price and physical regressions,
following the two time conventions used in the project. 
Before running this script, make sure you update the daily raw data excel file 
to cover the desired time period.

Outputs in "./bloomberg/":												   
1. bloomberg_prices_indep.frithu.xls - used for creating RHS in price regs.  	   
2. bloomberg_physical_indep.tue.xls - used for creating RHS in physical regs.
3. bloomberg_prices_dep.fri.xls - used for creating LHS in price regs.
4. bloomberg_prices_dep.mon.xls - used for creating LHS in price regs. (RDS)						

*/
* ---------------------------------------------------------------------------- *
***************************************************
****  (1) INDEPENDENT   ***************************
*** Prices regressions: Friday/Thursday Convention
*** Physical regressions: Tuesday Convention
***************************************************


*******************************************************************
* Friday/Thursday timing convention- to be used as RHS in prices regression
*******************************************************************

* We fill in missing Friday/Thursday observation with Thursday/Wednesday observation first. 
* If still missing, we fill in with Wednesday/Tuesday observation.
* For rdsaus, we are safe to use Friday measure. rdsa market ends Friday morning in EST due to timezone difference.
* For xomus and bpus(ADR), we use Thursday measure because US Stock market ends later than our 2:30pm cutoff for Futures market close. 
* For vol_30d we are not sure how this is calculated: whether it uses settlement price or last trade price. So we just use Thursday to be safe. 

* Format Friday data
* rdsaus
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

keep date rdsana eurusd

gen dow = dow(date)

gen rdsaus = rdsana * eurusd 

drop if dow != 5
drop dow rdsana eurusd
drop in 1
drop if _n > 1147

gen week = _n

rename rdsaus rdsaus_fri

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Thursday observations if Friday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

keep date rdsana eurusd
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 4
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename date date_thu
rename rdsaus rdsaus_thu

gen week = _n

* Merge Friday and Thursday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Wednesday observations if Friday & Thursday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

keep date rdsana eurusd
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 3
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename rdsaus rdsaus_wed
rename date date_wed

gen week = _n

* Merge Friday, Thursday, and Wednesday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

gen rdsaus = rdsaus_fri
replace rdsaus = rdsaus_thu if rdsaus == .
replace rdsaus = rdsaus_wed if rdsaus == .

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp_fri.dta", replace

* Format Thursday data (bpus, xomus, vol_30d)

import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 4
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename date date_thu
rename xomus xomus_thu
rename bpus bpus_thu
rename vol_30d vol_30d_thu

gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Wednesday observations if Thursday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 3
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename date date_wed
rename xomus xomus_wed
rename bpus bpus_wed
rename vol_30d vol_30d_wed

gen week = _n

* Merge Thursday and Wednesday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Tuesday observations if Thursday & Wednesday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 2
drop dow rdsana eurusd 

drop if _n > 1147

rename date date_tue
rename xomus xomus_tue
rename bpus bpus_tue
rename vol_30d vol_30d_tue

gen week = _n

* Merge Thursday, Wednesday, Tuesday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

gen xomus = xomus_thu
replace xomus = xomus_wed if xomus == .
replace xomus = xomus_tue if xomus == .

gen bpus = bpus_thu
replace bpus = bpus_wed if bpus == .
replace bpus = bpus_tue if bpus == .

gen vol_30d = vol_30d_thu
replace vol_30d = vol_30d_wed if vol_30d == .
replace vol_30d = vol_30d_tue if vol_30d == .

merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp_fri.dta", nogen

keep date xomus bpus vol_30d rdsaus  


export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_indep.frithu.xls", firstrow(variables) replace

*******************************************************************
* Tuesday timing convention - to be used as RHS of phyiscal regression
*******************************************************************
* Format Tuesday data
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 
drop if dow != 2

drop if _n > 1147

keep date vol_30d xomus bpus rdsaus

rename vol_30d vol_30d_tue
rename xomus xomus_tue
rename bpus bpus_tue
rename rdsaus rdsaus_tue

gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Monday observations if Tuesday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 1

drop if _n > 1147

keep date vol_30d xomus bpus rdsaus

rename date date_mon
rename vol_30d vol_30d_mon
rename xomus xomus_mon
rename bpus bpus_mon
rename rdsaus rdsaus_mon

gen week = _n

* Merge Tuesday and Monday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use last Friday's observations if Tuesday & Monday are missing

import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 
drop if dow != 5

drop if _n > 1147

keep date vol_30d xomus bpus rdsaus

rename date date_fri
rename vol_30d vol_30d_fri
rename xomus xomus_fri
rename bpus bpus_fri
rename rdsa rdsaus_fri

gen week = _n

* Merge Tuesday and Monday and Friday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

gen vol_30d = vol_30d_tue
replace vol_30d = vol_30d_mon if vol_30d == .
replace vol_30d = vol_30d_fri if vol_30d == .

gen xomus = xomus_tue
replace xomus = xomus_mon if xomus == .
replace xomus = xomus_fri if xomus == .

gen bpus = bpus_tue
replace bpus = bpus_mon if bpus == .
replace bpus = bpus_fri if bpus == .

gen rdsaus = rdsaus_tue
replace rdsaus = rdsaus_mon if rdsaus == .
replace rdsaus = rdsaus_fri if rdsaus == .


keep date vol_30d bpus xomus rdsaus

export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_physical_indep.tue.xls", firstrow(variables) replace


******************************************************
*****(2) Dependent*****Just for prices regression*****
******************************************************

*******************************************************************
* Friday timing convention
*******************************************************************

* We fill in missing Friday observation with Monday observation first. 
* If still missing, we fill in with Tuesday observation.

* Format Friday data
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 5
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

gen week = _n

rename xomus xomus_fri
rename bpus bpus_fri
rename vol_30d vol_30d_fri

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Monday observations if Friday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 1
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename date date_mon
rename xomus xomus_mon
rename bpus bpus_mon
rename vol_30d vol_30d_mon

gen week = _n

* Merge Friday and Monday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Tuesday observations if Friday & Monday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)

drop if dow != 2
drop dow rdsana eurusd 
drop in 1
drop if _n > 1147

rename date date_tue
rename xomus xomus_tue
rename bpus bpus_tue
rename vol_30d vol_30d_tue

gen week = _n

* Merge Friday, Monday, and Tuesday data
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

gen xomus = xomus_fri
replace xomus = xomus_mon if xomus == .
replace xomus = xomus_tue if xomus == .

gen bpus = bpus_fri
replace bpus = bpus_mon if bpus == .
replace bpus = bpus_tue if bpus == .

gen vol_30d = vol_30d_fri 
replace vol_30d = vol_30d_mon if vol_30d == .
replace vol_30d = vol_30d_tue if vol_30d == . 

keep date xomus bpus vol_30d

export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_dep.fri.xls", firstrow(variables) replace

*******************************************************************************************
* Monday RDS for Friday timing convention- We use monday measure due to timezone difference
*******************************************************************************************
* Format Monday data
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 1
keep date rdsaus
drop in 1
drop if _n > 1147

rename rdsaus rdsaus_mon
gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Tuesday observations if Monday is missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 2
keep date rdsaus
drop in 1
drop if _n > 1147

rename date date_tue
rename rdsaus rdsaus_tue
gen week = _n

* Merge 
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", replace 

* Use Wednesday observations if Monday and Tuesday are missing
import excel "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_raw.d.xlsx", sheet("ToStata") firstrow clear

drop G H
gen dow = dow(date)
gen rdsaus = rdsana * eurusd 

drop if dow != 3
keep date rdsaus
drop in 1
drop in 1 /*not a typo, we need it twice for wednesday to make 4/15 the first week*/
drop if _n > 1147

rename date date_wed
rename rdsaus rdsaus_wed
gen week = _n

* Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/temp.dta", nogen

gen rdsaus = rdsaus_mon
replace rdsaus = rdsaus_tue if rdsaus == .
replace rdsaus = rdsaus_wed if rdsaus == .

keep date rdsaus  

export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/bloomberg/bloomberg_prices_dep.mon.xls", firstrow(variables) replace

