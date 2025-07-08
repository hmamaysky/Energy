/*
This .do file based on previous final version of dataset to create longer horizons of variables.
Longer horizons: 24 weeks/ 52 weeks

Outputs:
1. transformed_data_prices_v19.4_mod_full_sample.dta
2. transformed_data_prices_v19.4_mod_restricted_sample.dta
3. transformed_data_physical_v19.4_mod_full_sample.dta
4. transformed_data_physical_v19.4_mod_restricted_sample.dta
*/

********************************************************************************
* Data Constructions: prices regressions (full-sample dataset)
********************************************************************************
use "C:\Users\13917\Desktop\CBS\Energy-RA\Fwd_ In sample results\transformed_data_prices_v19.3_mod_full_sample.dta", clear

*** Tells Stata the data are a time-series and seperate by a week
tsset date_Fri, delta(7)

*** Generate longer horizons for FutRet
* This is how FutRet is originally computed:
* gen FutRet_dep = (clc01 - L.clc01)/L.clc01
* replace FutRet_dep = (clc02 - L.clc02)/L.clc02 if newContract == 1 
* gen FutRet_t4 = (1 + FutRet_dep) * (1 + L.FutRet_dep) * (1 + L2.FutRet_dep) * (1 + L3.FutRet_dep) * 100

* Test if the implementation is correct by constructing 8-week
* gen FutRet_t8_test = (FutRet_t4_Fri * L4.FutRet_t4_Fri) / 100
* gen diff = abs(FutRet_t8_test - FutRet_t8_Fri)

* 24-week horizons						   
gen double FutRet_t24_Fri = (FutRet_t8_Fri * L8.FutRet_t8_Fri * L16.FutRet_t8_Fri) / 100^2

* 52-week horizons
gen double FutRet_t52_Fri = (FutRet_t4_Fri *  L4.FutRet_t4_Fri *  L8.FutRet_t4_Fri ///
                                   * L12.FutRet_t4_Fri * L16.FutRet_t4_Fri ///
                                   * L20.FutRet_t4_Fri * L24.FutRet_t4_Fri ///
                                   * L28.FutRet_t4_Fri * L32.FutRet_t4_Fri ///
                                   * L36.FutRet_t4_Fri * L40.FutRet_t4_Fri ///
                                   * L44.FutRet_t4_Fri * L48.FutRet_t4_Fri) / 100^12
								  
*** Generate longer horizons for DSpot, xomRet, bpRet, rdsaRet, DOilVol
* This is how DSpot is originally computed:
* gen DSpot_t4_Fri = log(price/L4.price) * 100

* This is how DOilVol is originally computed:
* DOilVol_t4_Fri = vol_30d - L4.vol_30d

* That's why we are using addition here

* Test if the implementation is correct by constructing 8-week
* gen DSpot_t8_test = DSpot_t4_Fri + L4.DSpot_t4_Fri
* gen diff = abs(DSpot_t8_test - DSpot_t8_Fri)
* gen DOilVol_t8_test = DOilVol_t4_Fri + L4.DOilVol_t4_Fri
* gen diff = abs(DOilVol_t8_test - DOilVol_t8_Fri)

* 24-week horizons
gen DSpot_t24_Fri = DSpot_t8_Fri + L8.DSpot_t8_Fri + L16.DSpot_t8_Fri
gen xomRet_t24_Fri = xomRet_t8_Fri + L8.xomRet_t8_Fri + L16.xomRet_t8_Fri
gen bpRet_t24_Fri = bpRet_t8_Fri + L8.bpRet_t8_Fri + L16.bpRet_t8_Fri
gen rdsaRet_t24_Mon = rdsaRet_t8_Mon + L8.rdsaRet_t8_Mon + L16.rdsaRet_t8_Mon
gen DOilVol_t24_Fri = DOilVol_t8_Fri + L8.DOilVol_t8_Fri + L16.DOilVol_t8_Fri

* 52-week horizons
foreach v in DSpot xomRet bpRet DOilVol{
    gen `v'_t52_Fri = ///
          `v'_t4_Fri  + L4.`v'_t4_Fri  + L8.`v'_t4_Fri  + L12.`v'_t4_Fri ///
        + L16.`v'_t4_Fri + L20.`v'_t4_Fri + L24.`v'_t4_Fri + L28.`v'_t4_Fri ///
        + L32.`v'_t4_Fri + L36.`v'_t4_Fri + L40.`v'_t4_Fri + L44.`v'_t4_Fri ///
        + L48.`v'_t4_Fri
}

gen rdsaRet_t52_Mon = ///
      rdsaRet_t4_Mon  + L4.rdsaRet_t4_Mon  + L8.rdsaRet_t4_Mon  + L12.rdsaRet_t4_Mon ///
    + L16.rdsaRet_t4_Mon + L20.rdsaRet_t4_Mon + L24.rdsaRet_t4_Mon + L28.rdsaRet_t4_Mon ///
    + L32.rdsaRet_t4_Mon + L36.rdsaRet_t4_Mon + L40.rdsaRet_t4_Mon + L44.rdsaRet_t4_Mon ///
    + L48.rdsaRet_t4_Mon

*** Generate longer horizons for WIPI
* 24 week horizons
gen WIPI_24wk_monthly = L16.WIPI_8wk_monthly
replace WIPI_24wk_monthly = WIPI_24wk_monthly[17] if _n == 16

* 52 week horizons
gen WIPI_52wk_monthly = L44.WIPI_8wk_monthly
replace WIPI_52wk_monthly = WIPI_52wk_monthly[45] if _n == 44

save "C:\Users\13917\Desktop\CBS\Energy-RA\Longer Horizon\transformed_data_prices_v19.4_mod_full_sample.dta", replace

********************************************************************************
* Data Constructions: prices regressions (restricted-sample dataset)
********************************************************************************
use "C:\Users\13917\Desktop\CBS\Energy-RA\Fwd_ In sample results\transformed_data_prices_v19.3_mod_restricted_sample.dta", clear

*** Tells Stata the data are a time-series and seperate by a week
tsset date_Fri, delta(7)

*** Generate longer horizons for FutRet
* 24-week horizons						   
gen double FutRet_t24_Fri = (FutRet_t8_Fri * L8.FutRet_t8_Fri * L16.FutRet_t8_Fri) / 100^2

* 52-week horizons
gen double FutRet_t52_Fri = (FutRet_t4_Fri *  L4.FutRet_t4_Fri *  L8.FutRet_t4_Fri ///
                                   * L12.FutRet_t4_Fri * L16.FutRet_t4_Fri ///
                                   * L20.FutRet_t4_Fri * L24.FutRet_t4_Fri ///
                                   * L28.FutRet_t4_Fri * L32.FutRet_t4_Fri ///
                                   * L36.FutRet_t4_Fri * L40.FutRet_t4_Fri ///
                                   * L44.FutRet_t4_Fri * L48.FutRet_t4_Fri) / 100^12
								  
*** Generate longer horizons for DSpot, xomRet, bpRet, rdsaRet, DOilVol
* 24-week horizons
gen DSpot_t24_Fri = DSpot_t8_Fri + L8.DSpot_t8_Fri + L16.DSpot_t8_Fri
gen xomRet_t24_Fri = xomRet_t8_Fri + L8.xomRet_t8_Fri + L16.xomRet_t8_Fri
gen bpRet_t24_Fri = bpRet_t8_Fri + L8.bpRet_t8_Fri + L16.bpRet_t8_Fri
gen rdsaRet_t24_Mon = rdsaRet_t8_Mon + L8.rdsaRet_t8_Mon + L16.rdsaRet_t8_Mon
gen DOilVol_t24_Fri = DOilVol_t8_Fri + L8.DOilVol_t8_Fri + L16.DOilVol_t8_Fri

* 52-week horizons
foreach v in DSpot xomRet bpRet DOilVol{
    gen `v'_t52_Fri = ///
          `v'_t4_Fri  + L4.`v'_t4_Fri  + L8.`v'_t4_Fri  + L12.`v'_t4_Fri ///
        + L16.`v'_t4_Fri + L20.`v'_t4_Fri + L24.`v'_t4_Fri + L28.`v'_t4_Fri ///
        + L32.`v'_t4_Fri + L36.`v'_t4_Fri + L40.`v'_t4_Fri + L44.`v'_t4_Fri ///
        + L48.`v'_t4_Fri
}

gen rdsaRet_t52_Mon = ///
      rdsaRet_t4_Mon  + L4.rdsaRet_t4_Mon  + L8.rdsaRet_t4_Mon  + L12.rdsaRet_t4_Mon ///
    + L16.rdsaRet_t4_Mon + L20.rdsaRet_t4_Mon + L24.rdsaRet_t4_Mon + L28.rdsaRet_t4_Mon ///
    + L32.rdsaRet_t4_Mon + L36.rdsaRet_t4_Mon + L40.rdsaRet_t4_Mon + L44.rdsaRet_t4_Mon ///
    + L48.rdsaRet_t4_Mon

*** Generate longer horizons for WIPI
* 24 week horizons
gen WIPI_24wk_monthly = L16.WIPI_8wk_monthly
replace WIPI_24wk_monthly = WIPI_24wk_monthly[17] if _n == 16

* 52 week horizons
gen WIPI_52wk_monthly = L44.WIPI_8wk_monthly
replace WIPI_52wk_monthly = WIPI_52wk_monthly[45] if _n == 44

save "C:\Users\13917\Desktop\CBS\Energy-RA\Longer Horizon\transformed_data_prices_v19.4_mod_restricted_sample.dta", replace

********************************************************************************
* Data Constructions: physical regressions (full-sample dataset)
********************************************************************************
use "C:\Users\13917\Desktop\CBS\Energy-RA\Fwd_ In sample results\transformed_data_physical_v19.3_mod_full_sample.dta", clear

*** Tells Stata the data are a time-series and seperate by a week
tsset date_Wed, delta(7)

*** Generate longer horizons for DProd
* Test if the implementation is correct by constructing 8-week
* local ln2_100 = 100*ln(2)
* gen DProd_t8_test = 100*ln(exp(DProd_t4_Wed/100)+1)  +  L4.DProd_t4_Wed  -  `ln2_100'
* gen diff = abs(DProd_t8_test - DProd_t8_Wed)

* 24 week horizons
tempvar lnT lnCum
gen double r_0 = DProd_t4_Wed/100              // r_0

forvalues k = 1/5 {                           // r_1 … r_5
    gen double r_`k' = L`=4*`k''.DProd_t4_Wed/100
}

gen double `lnCum' = r_5                       // c_5
gen double `lnT'   = exp(`lnCum')             // T = e^{c_5}

forvalues k = 4(-1)0 {                        // c_4 … c_0
    replace `lnCum' = r_`k' + `lnCum'
    replace `lnT'   = `lnT' + exp(`lnCum')
}

gen DProd_t24_Wed = (ln(`lnT') - ln(6)) * 100

* 52 week horizons
forvalues k = 6/12 {
    gen double r_`k' = L`=4*`k''.DProd_t4_Wed/100   // r_6 … r_12
}

replace `lnCum' = r_12
replace `lnT'   = exp(`lnCum')

forvalues k = 11(-1)0 {
    replace `lnCum' = r_`k' + `lnCum'
    replace `lnT'   = `lnT' + exp(`lnCum')
}

gen DProd_t52_Wed = (ln(`lnT') - ln(13)) * 100

* Clean-up
ds r_*
drop `r(varlist)'
drop `lnCum' `lnT'

*** Generate longer horizons for DInv
* Test if the implementation is correct by constructing 8-week
* gen DInv_t8_test = DInv_t4_Wed + L4.DInv_t4_Wed
* gen diff = abs(DInv_t8_test - DInv_t8_Wed)

* 24 week horizons
gen DInv_t24_Wed = DInv_t8_Wed + L8.DInv_t8_Wed + L16.DInv_t8_Wed

* 52 week horizons
gen DInv_t52_Wed = ///
      DInv_t4_Wed  + L4.DInv_t4_Wed  + L8.DInv_t4_Wed  + L12.DInv_t4_Wed ///
    + L16.DInv_t4_Wed + L20.DInv_t4_Wed + L24.DInv_t4_Wed + L28.DInv_t4_Wed ///
    + L32.DInv_t4_Wed + L36.DInv_t4_Wed + L40.DInv_t4_Wed + L44.DInv_t4_Wed ///
    + L48.DInv_t4_Wed
	
*** Generate longer horizons for WIPI
* 24 week horizons
gen WIPI_24wk_monthly = L16.WIPI_8wk_monthly
replace WIPI_24wk_monthly = WIPI_24wk_monthly[17] if _n == 16

* 52 week horizons
gen WIPI_52wk_monthly = L44.WIPI_8wk_monthly
replace WIPI_52wk_monthly = WIPI_52wk_monthly[45] if _n == 44
save "C:\Users\13917\Desktop\CBS\Energy-RA\Longer Horizon\transformed_data_physical_v19.4_mod_full_sample.dta", replace

********************************************************************************
* Data Constructions: physical regressions (restricted-sample dataset)
********************************************************************************
use "C:\Users\13917\Desktop\CBS\Energy-RA\Fwd_ In sample results\transformed_data_physical_v19.3_mod_restricted_sample.dta", clear

*** Tells Stata the data are a time-series and seperate by a week
tsset date_Wed, delta(7)

*** Generate longer horizons for DProd
* 24 week horizons
tempvar lnT lnCum
gen double r_0 = DProd_t4_Wed/100              // r_0

forvalues k = 1/5 {                           // r_1 … r_5
    gen double r_`k' = L`=4*`k''.DProd_t4_Wed/100
}

gen double `lnCum' = r_5                       // c_5
gen double `lnT'   = exp(`lnCum')             // T = e^{c_5}

forvalues k = 4(-1)0 {                        // c_4 … c_0
    replace `lnCum' = r_`k' + `lnCum'
    replace `lnT'   = `lnT' + exp(`lnCum')
}

gen DProd_t24_Wed = (ln(`lnT') - ln(6)) * 100

* 52 week horizons
forvalues k = 6/12 {
    gen double r_`k' = L`=4*`k''.DProd_t4_Wed/100   // r_6 … r_12
}

replace `lnCum' = r_12
replace `lnT'   = exp(`lnCum')

forvalues k = 11(-1)0 {
    replace `lnCum' = r_`k' + `lnCum'
    replace `lnT'   = `lnT' + exp(`lnCum')
}

gen DProd_t52_Wed = (ln(`lnT') - ln(13)) * 100

* Clean-up
ds r_*
drop `r(varlist)'
drop `lnCum' `lnT'

*** Generate longer horizons for DInv
* 24 week horizons
gen DInv_t24_Wed = DInv_t8_Wed + L8.DInv_t8_Wed + L16.DInv_t8_Wed

* 52 week horizons
gen DInv_t52_Wed = ///
      DInv_t4_Wed  + L4.DInv_t4_Wed  + L8.DInv_t4_Wed  + L12.DInv_t4_Wed ///
    + L16.DInv_t4_Wed + L20.DInv_t4_Wed + L24.DInv_t4_Wed + L28.DInv_t4_Wed ///
    + L32.DInv_t4_Wed + L36.DInv_t4_Wed + L40.DInv_t4_Wed + L44.DInv_t4_Wed ///
    + L48.DInv_t4_Wed
	
*** Generate longer horizons for WIPI
* 24 week horizons
gen WIPI_24wk_monthly = L16.WIPI_8wk_monthly
replace WIPI_24wk_monthly = WIPI_24wk_monthly[17] if _n == 16

* 52 week horizons
gen WIPI_52wk_monthly = L44.WIPI_8wk_monthly
replace WIPI_52wk_monthly = WIPI_52wk_monthly[45] if _n == 44
save "C:\Users\13917\Desktop\CBS\Energy-RA\Longer Horizon\transformed_data_physical_v19.4_mod_restricted_sample.dta", replace