/* 
This .do file aggregates the text measures from daily to weekly frequency.    
It performs the aggregation twice, once assuming Friday is the last day of   
the week, and another assuming Tuesday is the last day of the week. The 	   
former is for regressions explaining price variables, and the latter for 	   
those explaining variation in production or inventories (which I refer to    
as physical variables). 													   
																			   
Outputs in "./text_measures/":												   
1. text_measures_prices.fri - used in price regs.; Friday end-day 		   		   
2. text_measures_physical.tue - used in physical regs.; Tuesday end-day	   
																			   
The aggregation is a simple average of the observations in the 5 days of 	   
the week, depending on what we assume the last day of the week is.	       

*/
* ---------------------------------------------------------------------------- *
* PART I: Data cleaning common to both aggregations							   *
* ---------------------------------------------------------------------------- *

*** Import raw daily text measures data

import delimited "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_raw.d.csv", varnames(1) clear 

*** Rename varirables and drop unused ones 
forvalues i = 1/7 {
	rename topic`i' ftopic`i'
	rename topicsentiment`i' stopic`i'
}

rename articlecount artcount
drop unclassified unclassifiedsentiment

*** Format date
tostring date, gen(Date)
gen DATE = date(Date, "YMD")
format DATE %td
order DATE, before(Date)

*** Generate daily id variable 
gen dailydate = date(Date, "YMD")
gen dow 	  = dow(dailydate)

*** Check we only have obs. for Mon.-Fri.
tab dow
drop if dow == 0
tab dow 

*** Save data used for aggregation 
save "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/ToAggregate.dta", replace 


* ---------------------------------------------------------------------------- *
* PART II: Text measures aggregated assuming a Friday end-day				   *
* ---------------------------------------------------------------------------- *
use "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/ToAggregate.dta", clear 

*** Create weekly id using Sunday's weekid for a given week
gen weekid = dailydate - dow

*** Check that dow = 0
gen check = dow(weekid)
tab check 

*** Drop check if tab shows check = 0 only
drop check 

*** Drop date vars as appropriate, need complete weeks
*** Check beginning and end of dow, along with dow
*** to determine if we have incomplete weeks, then drop
drop if _n == 3
drop if _n == 2
drop if _n == 1

drop if _n == 5738
drop if _n == 5737
drop if _n == 5736

*** Clean up data for aggregation
drop date Date dailydate dow 
rename DATE date 
order weekid, before(artcount)
sort date

*** Convert to weekly series by average over days in a given week 
collapse (mean) artcount entropy ftopic1 ftopic2 ftopic3 ftopic4 ftopic5 ///
				ftopic6 ftopic7 stopic1 stopic2 stopic3 stopic4 stopic5  ///
				stopic6 stopic7 (last) date , by(weekid)

drop weekid
order date, before(artcount)

*** Save weekly measures 
export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_prices.fri_v3.xlsx", firstrow(variables) missing(".") replace


* ---------------------------------------------------------------------------- *
* PART III: Text measures aggregated assuming a Tuesday end-day				   *
* ---------------------------------------------------------------------------- *
use "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/ToAggregate.dta", clear 

*** Create weekly id assuming Tues. is the last day of the week
*** To do so, I create an adjustment factor series that would match
*** each day's id to earliest Tuesday's id (earlist ahead not pevious)
*** once added to dailydate
set obs 5741
egen adj = fill(6 5 4 1 0 6 5 4 1 0)

gen weekid = dailydate + adj

*** Check that dow variables in the order 3 4 5 1 2 share the same weekid, 
*** corresponding to dow() = 2
gen check = dow(weekid)
tab check 

*** Drop check if tab shows check = 2 only
drop check 

*** Drop date vars as appropriate
drop in 5741

*** Clean up data for aggregation
drop date Date dailydate dow 
rename DATE date 
order weekid, before(artcount)
sort date

*** Convert to weekly series by average over days in a given week 
collapse (mean) artcount entropy ftopic1 ftopic2 ftopic3 ftopic4 ftopic5 ///
				ftopic6 ftopic7 stopic1 stopic2 stopic3 stopic4 stopic5  ///
				stopic6 stopic7 (last) date , by(weekid)

drop weekid
order date, before(artcount)

*** save weekly measures 
export excel using "/Users/Economist/Dropbox/Research/ncm_research/data/text_measures/text_measures_physical.tue_v3.xlsx", firstrow(variables) missing(".") replace
