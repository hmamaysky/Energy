
# Energy
##  Variable construction

The variable construction code is described here.

### Descriptions for Codes
- ```update_bloomberg.do``` reads in daily Bloomberg data saved in data/bloomberg/bloomberg raw.d.xlsx, and creates follwoing four output files with weekly data for price and physical regressions, following the two time conventions used in the project. Before running this script, make sure you update the daily raw data excel file to cover the desired time period.

	- Outputs in "./bloomberg/":  

  		1. bloomberg_prices_indep.frithu.xls - used for creating RHS in price regs.  	   
  		2. bloomberg_physical_indep.tue.xls - used for creating RHS in physical regs.
  		3. bloomberg_prices_dep.fri.xls - used for creating LHS in price regs.
  		4. bloomberg_prices_dep.mon.xls - used for creating LHS in price regs. (RDS)					 

- ```update_fame.do``` reads in daily fame data saved in data/fame/fmae_raw.d.csv, and coverts to weekly series for both time conventions.   

	- Outputs in "./fame/":	 
	
  		1. fame_prices.thu_v4.csv - used for creating RHS in price regs.
  		2. fame_physical.tue_v4.csv - used for creating RHS in physical regs.
  		3. fame_prices.fri_v4.csv - contains the following: 
	   	(1) Friday measures for constructing LHS variables in Friday timing
	   	(2) Weekly prod and inv variables to construct DProd and DInv as RHS of Prices regressions

- ```update_futures_rea.do``` computes oil futures prices returns for both timing conventions. Code includes details on how we fill in missing observations. This script also includes the calculation of the World Industrial Production Index variable we use in our models.   

	- Outputs in "./futures_rea/":	 		   
	
		1. futures_rea_prices_indep.dta - used for RHS in price regs.  	   
		2. futures_rea_physical_indep.dta - used for RHS in physical regs.
		3. futures_rea_prices_dep - used for LHS in price regs.

- ```update_text_measures.do``` aggregates the text measures from daily to weekly frequency. It performs the aggregation twice, once assuming Friday is the last day of the week, and another assuming Tuesday is the last day of the week. The former is for regressions explaining price variables, and the latter for those explaining variation in production or inventories (which we refer to as physical variables). 		 										 

	- Outputs in "./text_measures/":    
										   
  		1. text_measures_prices.fri - used in price regs.; Friday end-day 		   		   
  		2. text_measures_physical.tue - used in physical regs.; Tuesday end-day	  

- ```merge_data.do``` does the following: 
  1. creates a "week" variable to identify weekly observations since within the same week obs. fall on diff. days for diff. vars
  2. creates four merged datasets, two following the time conventions used in oil price returns, oil company stock returns, and oil volatility regressions, and two others for oil prod. and inv.
  3. The datasets are in ./clean_data/ labeled for prices and physical
  4. The data in the merged_data .dta files is the underlying data for all transformations, thus, only transformations need to change when considering different time horizon   
		 										 

	- Outputs in "./clean_data/":     
		
  		1. merged_data_prices_dep.dta
  		2. merged_data_prices_indep.dta
  		3. merged_data_physical_dep.dta
  		4. merged_data_physical_indep.data

- ```transform_data.do``` performs the following: 
  1. transforms underlying merged data as needed for regressions
  2. creates 2 versions of the transformed data labeled prices or physical depending on the regressions it's meant to serve      

	- Outputs in "./clean_data/":	  
	
		1. transformed_data_prices_v19.dta
		2. transformed_data_physical_v19.dta

** Note that these are not the final versions of input data for the regressions because we add additional variables in subsequent .do files.

- ```va_weekly.do``` creates the following weekly variables.
  1. OpenInt
  2. liquidity 
  3. HedgPres   

	- Output in "./variable augmentation/":	    
	
	  	1. oi_tue.dta - used for RHS in physical regs
  		2. oi_thu.dta - used for RHS in price regs
  		3. liq_tue.dta - used for RHS in physical regs
  		4. liq_thu.dta - used for RHS in price regs
  		5. hp_fri.dta - used for RHS in price & physical regs

- ```va_to_python.do``` creates data that is exported to python for constructing Betas(InflaBeta, DolBeta)   

	- Output in "./variable augmentation/":	    
	
	  	1. cpidxyfutret_v2.dta 				 

Next, run python code ```va_beta.py``` to create InflaBeta and DolBeta    

- ```va_beta.py``` creates the following variables: 
  1. cpiyr_betas - renamed as InflaBeta in the final dataset
  2. dxy_betas - renamed as DolBeta in the final dataset       
   
	- Output in "./variable augmentation/":	    
	
	  	1. cpidxyfutret_v2.dta

- ```va_monthly.do``` creates the following monthly variables.
  1. BE/ME
  2. Mom
  3. BasMom 


	- Output in "./variable augmentation/":	     
	
  		1. data_va_monthly_v2.dta   

- ```va_merge.do``` incorporates newly created variables in ./variable augmentation/ to the final datasets.     

	- Output in "./clean_data/":	   
	
  		1. transformed_data_physical_v19.dta
		2. transformed_data_prices_v19.dta

- ```add_pca.py``` creates the following variables.
  1. PCAsent
  2. PCAfreq
  3. PCAall       
  
  	- Output in "./forwardSelection/":	   
	
  		1. transformed_data_physical_pca_v19.dta
  		2. transformed_data_prices_pca_v19.dta

- ```merge_risk_premia_pca.do``` transforms the relevant risk premia series for each timing convention, and merges them with the datasets of all other variables. In addition, this script incorporates pca variables created in python into the final datasets.       


	- Output in "./forwardSelection/":	   

  		1. transformed_data_physical_v19.dta
  		2. transformed_data_prices_v19.dta


### Notes before Running Codes
- Please change the working directory and the saving directory in each file appropriately before running the program.


### Procedures
- Please run the codes in the following order.
```
1. update_bloomberg.do
```
```
2. update_fame.do
```
```
3. update_futures_rea.do
```
```
4. update_text_measures.do
```
```
5. merge_data.do
```
```
6. transform_data.do
```
```
7. va_weekly.do
```
```
8. va_beta.py
```
```
9. va_monthly.do
```
```
10. va_merge.do
```
```
11. add_pca.py
```
```
12. merge_risk_premia_pca.do
```
