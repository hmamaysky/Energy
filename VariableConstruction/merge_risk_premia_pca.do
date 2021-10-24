/*
This .do file transforms the relevant risk premia series for each timing convention, 
and merges them with the datasets of all other variables. 
In addition, this script incorporates pca variables created in python into the final datasets.

Output in "./forwardSelection/":	
1. transformed_data_physical_v19.dta
2. transformed_data_prices_v19.dta

*/

******************************************************************************
* Risk premia following Thursday EOP for Fride EOP timing convention
******************************************************************************
********** Risk premia 
* Thursday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 4
drop if date < date("19980409", "YMD")
drop if date > date("20200326", "YMD")

gen vix_thu = VIX_px_last
gen spx_thu = SPX_vol_30d 
gen ovx_thu = OVX_px_last
gen cl1_thu = WTI_vol_30d

rename date date_Thu
keep date vix_thu spx_thu ovx_thu cl1_thu

gen week = _n

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", replace 

* Wednesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 3
drop if date < date("19980408", "YMD")
drop if date > date("20200325", "YMD")

gen vix_wed = VIX_px_last
gen spx_wed = SPX_vol_30d 
gen ovx_wed = OVX_px_last
gen cl1_wed = WTI_vol_30d

rename date date_wed

keep date_wed vix_wed spx_wed ovx_wed cl1_wed

gen week = _n

* Merge 
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", replace 

* Tuesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 2
drop if date < date("19980407", "YMD")
drop if date > date("20200324", "YMD")

gen vix_tue = VIX_px_last
gen spx_tue = SPX_vol_30d 
gen ovx_tue = OVX_px_last
gen cl1_tue = WTI_vol_30d

rename date date_tue

keep date_tue vix_tue spx_tue ovx_tue cl1_tue

gen week = _n

* Merge 
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", nogen

gen vix = vix_thu
replace vix = vix_wed if vix == .
replace vix = vix_tue if vix == .

gen spx = spx_thu
replace spx = spx_wed if spx == .
replace spx = spx_tue if spx == .

gen ovx = ovx_thu
replace ovx = ovx_wed if ovx == .
replace ovx = ovx_tue if ovx == .

gen cl1 = cl1_thu
replace cl1 = cl1_wed if cl1 == .
replace cl1 = cl1_tue if cl1 == .

gen vix_diff_Thu = vix - spx 
gen ovx_diff_Thu = ovx - cl1

keep week vix_diff_Thu ovx_diff_Thu

***** Merge newly created independent variables with existing series 
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_v19.dta", nogen
order vix_diff ovx_diff, before(BEME_monthly)
save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", replace 
*********************
***** Format SDF data
*********************
*Thursday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Thu
gen dow = dow(date_Thu)
drop if dow != 4

keep date_Thu E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_thu

* Merge SDF and risk premia 
merge 1:1 date_Thu using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Thu

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", replace 

* Wednesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Wed
gen dow = dow(date_Wed)
drop if dow != 3

keep date_Wed E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_wed

* Merge SDF and risk premia 
merge 1:1 date_Wed using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Wed

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", replace 

* Tuesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Tue
gen dow = dow(date_Tue)
drop if dow != 2

keep date_Tue E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_tue

* Merge SDF and risk premia 
merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Tue


gen sdf_fullSample_Thu = sdf_fullSample_thu
replace sdf_fullSample_Thu = sdf_fullSample_wed if sdf_fullSample_Thu == .
replace sdf_fullSample_Thu = sdf_fullSample_tue if sdf_fullSample_Thu == .

drop sdf_fullSample_thu sdf_fullSample_wed sdf_fullSample_tue

format date_Fri-date_Mon %tdnn/dd/CCYY

order week, after(monthly)
order sdf_fullSample_Thu, after(ovx_diff_Thu)

drop monthly week

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", replace

*********************************************
***** Add PCA variables to the final dataset
*********************************************
* import the pca data file created from python code 'add_pca.py'. make sure you change the directories as necessary.
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_prices_pca_v19.dta", clear

* Change date format to be compatible
gen temp = date(date_Fri, "YMD")
format %tdnn/dd/CCYY temp
drop date_Tue - date_Mon
rename temp date_Fri
order date_Fri
save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_pca_v19.dta", replace

* incorporate pca into the dataset

use "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", clear

merge 1:1 date_Fri using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_pca_v19.dta"
drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_prices_v19.dta", replace


******************************************************************************
* Risk premia following Tuesday EOP timing convention
******************************************************************************
******* Risk premia 
* Tuesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 2
drop if date < date("19980407", "YMD")
drop if date > date("20200324", "YMD")

gen vix_tue = VIX_px_last
gen spx_tue = SPX_vol_30d 
gen ovx_tue = OVX_px_last
gen cl1_tue = WTI_vol_30d

keep date vix_tue spx_tue ovx_tue cl1_tue

gen week = _n
rename date date_Tue

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", replace 

* Monday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 1
drop if date < date("19980406", "YMD")
drop if date > date("20200323", "YMD")

gen vix_mon = VIX_px_last
gen spx_mon = SPX_vol_30d 
gen ovx_mon = OVX_px_last
gen cl1_mon = WTI_vol_30d

keep date vix_mon spx_mon ovx_mon cl1_mon

gen week = _n
rename date date_mon

* Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", nogen

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", replace 

* Friday(last week)
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 5
drop if date < date("19980403", "YMD")
drop if date > date("20200320", "YMD")

gen vix_fri = VIX_px_last
gen spx_fri = SPX_vol_30d 
gen ovx_fri = OVX_px_last
gen cl1_fri = WTI_vol_30d

keep date vix_fri spx_fri ovx_fri cl1_fri

gen week = _n
rename date date_fri

* Merge
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/temp.dta", nogen

gen vix = vix_tue
replace vix = vix_mon if vix == .
replace vix = vix_fri if vix == .

gen spx = spx_tue
replace spx = spx_mon if spx == .
replace spx = spx_fri if spx == .

gen ovx = ovx_tue
replace ovx = ovx_mon if ovx == .
replace ovx = ovx_fri if ovx == .

gen cl1 = cl1_tue
replace cl1 = cl1_mon if cl1 == .
replace cl1 = cl1_fri if cl1 == .

gen vix_diff_Tue = vix - spx 
gen ovx_diff_Tue = ovx - cl1

keep date_Tue week vix_diff_Tue ovx_diff_Tue

******* Merge newly created independent variables with existing series 
merge 1:1 week using "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_v19.dta", nogen

order vix_diff_Tue ovx_diff_Tue, before(BEME_monthly)

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace 

******* Format SDF data
* Tuesday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Tue
gen dow = dow(date_Tue)
drop if dow != 2

keep date_Tue E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_tue

* Merge SDF and risk premia 
merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Tue

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace

* Monday
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Mon
gen dow = dow(date_Mon)
drop if dow != 1

keep date_Mon E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_mon

gen date_Tue = date_Mon+1
format date_Tue %tdnn/dd/CCYY

* Merge SDF and risk premia 
merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Tue

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace

* Friday(Prior)
import excel "/Users/Economist/Dropbox/Research/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_Fri
gen dow = dow(date_Fri)
drop if dow != 5

keep date_Fri E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_fri

gen date_Tue = date_Fri +4
format date_Tue %tdnn/dd/CCYY

* Merge SDF and risk premia 
merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta"

drop if _merge == 1
drop _merge 
sort date_Tue

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace


gen sdf_fullSample_Tue = sdf_fullSample_tue
replace sdf_fullSample_Tue = sdf_fullSample_mon if sdf_fullSample_Tue == .
replace sdf_fullSample_Tue = sdf_fullSample_fri if sdf_fullSample_Tue == .

drop sdf_fullSample_tue sdf_fullSample_mon sdf_fullSample_fri date_Fri date_Mon

order date_Tue date_Wed monthly week
order sdf_fullSample_Tue, after(ovx_diff_Tue)
order tnote_10y_Tue, after(VIX_Tue)
format date_Wed %tdnn/dd/CCYY

drop monthly week

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace


*********************************************
***** Add PCA variables to the final dataset
*********************************************
* import the pca data file created from python code 'add_pca.py'. make sure you change the directories as necessary.
use "/Users/Economist/Dropbox/Research/ncm_research/data/clean_data/transformed_data_physical_pca_v19.dta", clear

* Change date format to be compatible
gen temp = date(date_Tue, "YMD")
format %tdnn/dd/CCYY temp
drop date_Tue - date_Wed
rename temp date_Tue
order date_Tue

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_pca_v19.dta", replace

* incorporate pca into the dataset

use "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", clear

merge 1:1 date_Tue using "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_pca_v19.dta"
drop _merge

save "/Users/Economist/Dropbox/Research/ncm_research/forwardSelection/transformed_data_physical_v19.dta", replace

