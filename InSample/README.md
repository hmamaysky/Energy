# Energy
## In Sample Analysis

Contains codes for In Sample analysis for the Energy Project

### Descriptions for Codes and Data
- In ForwardSelection folder: 
  - stepwiseForwardSelection_detrended.R performs stepwise forward selection to choose seven regressors from all the candidates.
- In Monte Carlo folder:
  - Monte_Carlo_forwardSelection.R performs monte carlo simulation to control for small-sample biases in R-squareds and t-statistics.
- In Tables folder: 
  - Table_V.R: Table V- "Measuring Instability of In-Sample Forward Selection Results"
  - Table_A_IV.R: Table A.IV- "Count of Selected and Statistically Significant Variables of Forward Selection Model"
  - Table_A_V.R: Table V- "F-test for the Eight-Week Stepwise Forward Selection"
  - Table_A_VII.R: Table A.VII- "Subperiod stepwise forward selection at the eight-week horizon"
 
### Notes for Further Processing of Results
- Results for Forward Selection Model is self explanatory. Each specification has a word file for regression stats and the order of selected variables is the same as the row orders. The results are used for creating Table A.VI. 
- Results for Monte Carlo analysis are well formatted. We can combine them with Forward Selection results to create Table IV.
- Figure 3 and Figure 4 are created in Monte_Carlo_forwardSelection.R.
- R codes in Tables folder are used to create some tables in the paper. (Table V, Table A.IV, Table A.V, Table A.VII) 

### Notes before running Analysis
- Please change the working directory and the saving directory properly in each file before running the program.
- Run stepwiseForawrdSelection_detrended.R before running Table_V.R
- Run Monte_Carlo_forwardSelection.R before running Table_A_IV.R

### Procedures
1. In Sample Forward Selection Model
```
stepwiseForwardSelection_detrended.R
```

2. Monte Carlo Simulation
```
Monte_Carlo_forwardSelection.R
```
3. Creating Some Tables
```
Table_V.R
Table_A_IV.R 
Table_A_V.R
Table_A_VII.R
```
