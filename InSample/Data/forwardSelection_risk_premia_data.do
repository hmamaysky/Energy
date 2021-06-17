******************************************************************************
* Risk premia following Thursday EOP for Fride EOP timing convention
* Last Version: 2020-11-22
******************************************************************************
********** Risk premia 
* Thursday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

gen dow = dow(date)
drop if dow != 4
drop if date < date("19980409", "YMD")
drop if date > date("20200326", "YMD")

gen vix_thu = VIX_px_last
gen spx_thu = SPX_vol_30d 
gen ovx_thu = OVX_px_last
gen cl1_thu = WTI_vol_30d

keep date vix_thu spx_thu ovx_thu cl1_thu

gen week = _n

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/temp.dta", replace 

* Wednesday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

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
merge 1:1 week using "/Users/billwu/Dropbox/ncm_research/forwardSelection/temp.dta", nogen

gen vix = vix_thu
replace vix = vix_wed if vix == .

gen spx = spx_thu
replace spx = spx_wed if spx == .

gen ovx = ovx_thu
replace ovx = ovx_wed if ovx == .

gen cl1 = cl1_thu
replace cl1 = cl1_wed if cl1 == .

gen vix_spx = vix - spx 
gen ovx_cl1 = ovx - cl1

rename date date_thurs

keep date_thurs week vix_spx ovx_cl1

***** Merge newly created independent variables with existing series 
merge m:1 week using "/Users/billwu/Dropbox/ncm_research/data/clean_data/transformed_data_prices_v14.dta", nogen keep(match)

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_prices_v14.dta", replace 

***** Format SDF data
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_thurs
gen dow = dow(date)
drop if dow != 4

keep date_thurs E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_thu

* Merge SDF and risk premia 
merge m:1 date_thurs using "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_prices_v14.dta"

drop if _merge == 1
drop _merge 
sort date

gen date_wed = date_thurs - 1
format date_wed %tddd-Mon-YY

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_prices_v14.dta", replace 

* Wednesday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_wed
gen dow = dow(date)
drop if dow != 3

keep date_wed E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_wed

* Merge SDF and risk premia 
merge m:1 date_wed using "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_prices_v14.dta"

drop if _merge == 1
drop _merge 
sort date

gen sdf_fullSample = sdf_fullSample_thu
replace sdf_fullSample = sdf_fullSample_wed if sdf_fullSample == .

drop sdf_fullSample_thu sdf_fullSample_wed

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_prices_v14.dta", replace

******************************************************************************
* Risk premia following Tuesday EOP timing convention
******************************************************************************
******* Risk premia 
* Tuesday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

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
rename date date_tue

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/temp.dta", replace 

* Monday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/risk_premia_data.xlsx", sheet("Sheet1") firstrow clear

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
merge 1:1 week using "/Users/billwu/Dropbox/ncm_research/forwardSelection/temp.dta"

gen vix = vix_tue
replace vix = vix_mon if vix == .

gen spx = spx_tue
replace spx = spx_mon if spx == .

gen ovx = ovx_tue
replace ovx = ovx_mon if ovx == .

gen cl1 = cl1_tue
replace cl1 = cl1_mon if cl1 == .

gen vix_spx = vix - spx 
gen ovx_cl1 = ovx - cl1

keep date_tue week vix_spx ovx_cl1

******* Merge newly created independent variables with existing series 
merge m:1 week using "/Users/billwu/Dropbox/ncm_research/data/clean_data/transformed_data_physical_v14.dta", nogen keep(match)

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_physical_v14.dta", replace 

******* Format SDF data
* Tuesday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_tue
gen dow = dow(date)
drop if dow != 2

keep date_tue E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_tue

* Merge SDF and risk premia 
merge m:1 date_tue using "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_physical_v14.dta"

drop if _merge == 1
drop _merge 
sort date

gen date_mon = date_tue - 1
format date_mon %tddd-Mon-YY

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_physical_v14.dta", replace

* Monday
import excel "/Users/billwu/Dropbox/ncm_research/regressions/v11/SDF-full-sample-2020-04-13.xlsx", sheet("SDF-full-sample-2020-04-13") firstrow clear

rename A date_mon
gen dow = dow(date)
drop if dow != 1

keep date_mon E_wti_fut_tr
rename E_wti_fut_tr sdf_fullSample_mon

* Merge SDF and risk premia 
merge m:1 date_mon using "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_physical_v14.dta"

drop if _merge == 1
drop _merge 
sort date

gen sdf_fullSample = sdf_fullSample_tue
replace sdf_fullSample = sdf_fullSample_mon if sdf_fullSample == .

drop sdf_fullSample_tue sdf_fullSample_mon

save "/Users/billwu/Dropbox/ncm_research/forwardSelection/transformed_data_physical_v14.dta", replace
