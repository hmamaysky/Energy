# Energy
## In Sample Analysis

Contains codes for In Sample analysis for the Energy Project

### Descriptions for Codes and Data
- In ForwardSelection folder: 
  - ```stepwiseForwardSelection_detrended.R``` performs stepwise forward selection to choose seven regressors from all the candidates.
- In Monte Carlo folder:
  - ```Monte_Carlo_forwardSelection.R``` performs monte carlo simulation to control for small-sample biases in R-squareds and t-statistics.
  - Figure 3 and Figure 4 are created in Monte_Carlo_forwardSelection.R.

- In Tables folder: 
  - R codes in this folder are used to create some tables in the paper. 
    - ```Table_V.R```: Table V- "Measuring Instability of In-Sample Forward Selection Results"
    - ```Table_A_IV.R```: Table A.IV- "Count of Selected and Statistically Significant Variables of Forward Selection Model"
    - ```Table_A_V.R```: Table V- "F-test for the Eight-Week Stepwise Forward Selection"
    - ```Table_A_VII.R```: Table A.VII- "Subperiod stepwise forward selection at the eight-week horizon"
 
### Notes for Further Processing of Results
- Results of Forward Selection Model are self explanatory. Each specification has a word file for regression stats and the order of selected variables is the same as the row orders. The results are used for creating Table IV and Table A.VI. 
- We use results of ```Monte_Carlo_forwardSelection.R``` and ```stepwiseForwardSelection_detrended.R``` to create Table IV.


### Notes before running Analysis
- Please change the working directory and the saving directory properly in each file before running the program.
- Run ```stepwiseForawrdSelection_detrended.R``` before running ```Table_V.R```
- Run ```Monte_Carlo_forwardSelection.R``` before running ```Table_A_IV.R```

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
