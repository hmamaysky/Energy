# Energy
## In Sample Analysis

Contains codes for In Sample analysis for the Energy Project

### Descriptions for Codes and Data
- stepwiseForawrdSelection.R 
- Monte_Carlo_forwardSelection.R 
- Data Section contains codes for the construction for non-textual data, and the .py file to add PCA series to the textual data. (Need further work)

### Notes for Further Processing of Results
- Results for Forward Selection Model is self explanatory. Each specification has an excel for regression stats and the order of selected variables is the same as the row orders.
- Graphs and tables from Monte Carlo analysis are well formatted, we can integrate them manually (tiny work as most are done in the code) to create tables and graphs in the paper.
- To collect the standardized coefficients, refer to line 690-740 in Monte_Carlo_forwardSelection.R.

### Notes before running Analysis
- Please change the working directory and the saving directory in each file before running the program.
- Required Data can also be found on the research grid under folder shared/energy_drivers/2020-11-16/data

### Procedures
1. In Sample Forward Selection Model

2. Monte Carlo Simulation

