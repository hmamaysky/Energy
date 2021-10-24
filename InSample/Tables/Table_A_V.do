/*
This .do file creates Table A.V "F-test for the Eight-Week Stepwise Forward Selection".

*/

*****Price

use "C:\Users\Economist\Dropbox\Research\ncm_research\forwardSelection\transformed_data_prices_v19.dta"


tsset date_Fri, delta(7 days)

**** FutRet
eststo m1: newey FutRet_t8_Fri L8.FutRet_Fri L8.tnote_10y_Thu L8.BEME_monthly L8.Mom_monthly L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fBbl_Fri, lag(6) force

*text
test L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fBbl_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.tnote_10y_Thu L8.BEME_monthly L8.Mom_monthly
mat t2= r(F)
mat p2= r(p)
*all
test L8.tnote_10y_Thu L8.BEME_monthly L8.Mom_monthly L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fBbl_Fri
mat t3= r(F)
mat p3= r(p)

*record the results in a matrix for the purpose of creating a table in the end

mat m1= t1,t2,t3
mat colnames m1= text notext all
mat pval1=p1,p2,p3
mat colnames pval1= text notext all
estadd matrix m= m1
estadd matrix pval= pval1

**** DSpot 
eststo m2: newey DSpot_t8_Fri L8.DSpot_Fri L8.basis_Fri L8.Mom_monthly L8.VIX_Thu L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fRpc_Fri, lag(6) force

*text
test L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fRpc_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.basis_Fri L8.Mom_monthly L8.VIX_Thu
mat t2= r(F)
mat p2= r(p)
*all
test L8.basis_Fri L8.Mom_monthly L8.VIX_Thu L8.PCAall_Fri L8.entropy_Fri L8.sGom_Fri L8.fRpc_Fri
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m2= t1,t2,t3
mat colnames m2= text notext all
mat pval2=p1,p2,p3
mat colnames pval2= text notext all
estadd matrix m= m2
estadd matrix pval= pval2

**** DOilvol 
eststo m3: newey DOilVol_t8_Fri L8.DOilVol_Thu L8.DSpot_Fri L8.OilVol_Thu WIPI_8wk_monthly L8.VIX_Thu L8.entropy_Fri L8.fCo_Fri L8.fGom_Fri, lag(6) force

*text
test L8.entropy_Fri L8.fCo_Fri L8.fGom_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.DSpot_Fri L8.OilVol_Thu WIPI_8wk_monthly L8.VIX_Thu
mat t2= r(F)
mat p2= r(p)
*all
test L8.DSpot_Fri L8.OilVol_Thu WIPI_8wk_monthly L8.VIX_Thu L8.entropy_Fri L8.fCo_Fri L8.fGom_Fri
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m3= t1,t2,t3
mat colnames m3= text notext all
mat pval3=p1,p2,p3
mat colnames pval3= text notext all
estadd matrix m= m3
estadd matrix pval= pval3


**** xomRet
eststo m4: newey xomRet_t8_Fri L8.xomRet_Thu L8.DInv_Wed L8.BEME_monthly L8.entropy_Fri L8.HedgPres_Fri L8.sEpg_Fri L8.fBbl_Fri L8.sEp_Fri, lag(6) force

*text
test L8.entropy_Fri L8.sEpg_Fri L8.fBbl_Fri L8.sEp_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.DInv_Wed L8.HedgPres_Fri L8.BEME_monthly
mat t2= r(F)
mat p2= r(p)
*all
test L8.DInv_Wed L8.BEME_monthly L8.entropy_Fri L8.HedgPres_Fri L8.sEpg_Fri L8.fBbl_Fri L8.sEp_Fri
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m4= t1,t2,t3
mat colnames m4= text notext all
mat pval4=p1,p2,p3
mat colnames pval4= text notext all
estadd matrix m= m4
estadd matrix pval= pval4

**** bpRet 
eststo m5: newey bpRet_t8_Fri L8.bpRet_Thu L8.basis_Fri L8.tnote_10y_Thu L8.Mom_monthly L8.vix_diff_Thu L8.entropy_Fri L8.sGom_Fri L8.sEnv_Fri, lag(6) force

*text
test L8.entropy_Fri L8.sGom_Fri L8.sEnv_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.basis_Fri L8.tnote_10y_Thu L8.Mom_monthly L8.vix_diff_Thu
mat t2= r(F)
mat p2= r(p)
*all
test L8.basis_Fri L8.tnote_10y_Thu L8.Mom_monthly L8.vix_diff_Thu L8.entropy_Fri L8.sGom_Fri L8.sEnv_Fri
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m5= t1,t2,t3
mat colnames m5= text notext all
mat pval5=p1,p2,p3
mat colnames pval5= text notext all
estadd matrix m= m5
estadd matrix pval= pval5

**** rdsaRet 
eststo m6: newey rdsaRet_t8_Mon L8.rdsaRet_Fri L8.DInv_Wed L8.InflaBeta_monthly L8.Mom_monthly L8.entropy_Fri L8.fCo_Fri L8.sGom_Fri L8.sEnv_Fri, lag(6) force

*text
test L8.entropy_Fri L8.fCo_Fri L8.sGom_Fri L8.sEnv_Fri
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.DInv_Wed L8.InflaBeta_monthly L8.Mom_monthly
mat t2= r(F)
mat p2= r(p)
*all
test L8.DInv_Wed L8.InflaBeta_monthly L8.Mom_monthly L8.entropy_Fri L8.fCo_Fri L8.sGom_Fri L8.sEnv_Fri
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m6= t1,t2,t3
mat colnames m6= text notext all
mat pval6=p1,p2,p3
mat colnames pval6= text notext all
estadd matrix m= m6
estadd matrix pval= pval6


******Physical

use "C:\Users\Economist\Dropbox\Research\ncm_research\forwardSelection\transformed_data_physical_v19.dta"

tsset date_Tue, delta(7 days)

**** DInv 

eststo m7: newey DInv_t8_Wed L8.DInv_Wed L8.tnote_10y_Tue L8.HedgPres_Fri L8.entropy_Tue L8.sGom_Tue L8.fRpc_Tue L8.sEp_Tue L8.fEp_Tue , lag(6) force

*text
test L8.entropy_Tue L8.sGom_Tue L8.fRpc_Tue L8.sEp_Tue L8.fEp_Tue
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.tnote_10y_Tue L8.HedgPres_Fri
mat t2= r(F)
mat p2= r(p)
*all
test L8.tnote_10y_Tue L8.HedgPres_Fri L8.entropy_Tue L8.sGom_Tue L8.fRpc_Tue L8.sEp_Tue L8.fEp_Tue
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m7= t1,t2,t3
mat colnames m7= text notext all
mat pval7=p1,p2,p3
mat colnames pval7= text notext all
estadd matrix m= m7
estadd matrix pval= pval7

**** DProd

eststo m8: newey DProd_t8_Wed L8.DProd_Wed L8.DSpot_Tue L8.BEME_monthly L8.InflaBeta_monthly L8.VIX_Tue L8.fGom_Tue L8.fBbl_Tue L8.sRpc_Tue , lag(6) force

*text
test L8.fGom_Tue L8.fBbl_Tue L8.sRpc_Tue
mat t1= r(F)
mat p1= r(p)
*nontext
test L8.DSpot_Tue L8.BEME_monthly L8.InflaBeta_monthly L8.VIX_Tue
mat t2= r(F)
mat p2= r(p)
*all
test L8.DSpot_Tue L8.BEME_monthly L8.InflaBeta_monthly L8.VIX_Tue L8.fGom_Tue L8.fBbl_Tue L8.sRpc_Tue
mat t3= r(F)
mat p3= r(p)
*record the results in a matrix for the purpose of creating a table in the end
mat m8= t1,t2,t3
mat colnames m8= text notext all
mat pval8=p1,p2,p3
mat colnames pval8= text notext all
estadd matrix m= m8
estadd matrix pval= pval8

*create a single table with Fstatistics for the 8 dependent variables
esttab m* using C:\Users\Economist\Dropbox\Research\ncm_research\forwardSelection\20210821\Table_A_V.rtf, replace star(* 0.10 ** 0.05 *** 0.01) cells(m(star pvalue(pval) fmt(%9.2f))) ///
collab(none) mlab("FutRet" "DSpot" "DOilVol" "xomRet" "bpRet" "rdsaRet" "DInv" "DProd") nonumb noobs ///
coeflab(notext "non-text") note( "Note: *, **, and *** represent p<0.1, p<0.05, and p<0.01, respectively.")
