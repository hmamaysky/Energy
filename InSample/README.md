# Energy
## In Sample Analysis

Contains codes for In Sample analysis for the Energy Project

### Descriptions for Codes and Data
- In ForwardSelection folder: 
  - ```stepwiseForwardSelection_detrended.R``` performs stepwise forward selection to choose seven regressors from all the candidates.
- In Monte Carlo folder:
  - ```Monte_Carlo_forwardSelection.R``` performs monte carlo simulation to control for small-sample biases in R-squareds and t-statistics.
  - Figure 3 and Figure 4 are created in Monte_Carlo_forwardSelection.R.

- In Tables and Figures folder: 
  - R codes in this folder are used to create some tables in the paper. 
    - ```Table_V.R```: Table V- "Measuring Instability of In-Sample Forward Selection Results"
    - ```Table_A_IV.R```: Table A.IV- "Count of Selected and Statistically Significant Variables of Forward Selection Model"
    - ```Table_A_V.R```: Table V- "F-test for the Eight-Week Stepwise Forward Selection"
    - ```Table_A_VII.R```: Table A.VII- "Subperiod stepwise forward selection at the eight-week horizon"
    - ```Fig_1.py```: Figure 1- "Word cloud plots for energy topics" 
    - ```Fig_2.py```: Figure 2- "NLP measures over time"
    - ```Fig_A_1.py```: Figure A.1- "Clustered Correlation Plot of all Independent Variable Series"
 
 
### Notes for Further Processing of Results
- Results of Forward Selection Model are self explanatory. Each specification has a word file for regression stats and the order of selected variables is the same as the row orders. The results are used for creating Table IV and Table A.VI. 
- We utilize both the results of ```stepwiseForwardSelection_detrended.R``` and ```Monte_Carlo_forwardSelection.R``` to create Table IV. While information about the order of selected variables comes from the former, other results for creating Table IV such as coefficients, p-values, mean of simulated adjusted, CDF(%) come from the latter.

### Notes for running analysis with Michael Plante data
- files _opec, _mnsc, and _all run are all used to run this analysis for the appendix. _mnsc runs the same analysis as table 4 in the paper with the same data, _opec runs the same analysis as table 4 in the paper but without any of the text analysis variables and with Michael Plante's OPEC variable, and _all runs the same analysis as table 4 in the paper but with both the papers text analysis variables and with Michael Plante's OPEC variable.


### Notes before running Analysis
- Please change the working directory and the saving directory properly in each file before running the program.
- Run ```stepwiseForawrdSelection_detrended.R``` before running ```Table_V.R```
- Run ```Monte_Carlo_forwardSelection.R``` before running ```Table_A_IV.R```

### Procedures
1. In Sample Forward Selection Model
```
- stepwiseForwardSelection_detrended.R
```

2. Monte Carlo Simulation
```
- Monte_Carlo_forwardSelection.R
```
3. Creating Some Tables
```
- Table_V.R
- Table_A_IV.R 
- Table_A_V.R
- Table_A_VII.R
- Fig_1.py
- Fig_2.py
- Fig_A_1.py
```
