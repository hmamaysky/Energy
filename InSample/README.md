# Energy
## In Sample Analysis

Contains codes for In Sample analysis for the Energy Project

### Descriptions for Codes and Data
- In ForwardSelection folder: 
  - stepwiseForwardSelection_detrended.R 
- In Monte Carlo folder:
  - Monte_Carlo_forwardSelection.R 
- In Tables folder: 
  - Table_V.R
  - Table_A_IV.R
  - Table_A_V.R
  - Table_A_VII.R
 
### Notes for Further Processing of Results
- Results for Forward Selection Model is self explanatory. Each specification has an word file for regression stats and the order of selected variables is the same as the row orders. The results are used for creating Table A.VI. 
- Graphs and tables from Monte Carlo analysis are well formatted, we can integrate them manually (tiny work as most are done in the code) to create tables and graphs in the paper. (Table IV, Figure 3, Figure 4)
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
