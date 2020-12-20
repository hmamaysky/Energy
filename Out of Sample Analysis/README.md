# Energy
## Out of Sample Analysis

Contains codes for OOS analysis for the Energy Project

### Descriptions for Codes
#### Codes for Analysis
- OOSfuncs.py includes all the functions needed for the OOS analysis. 
- Forward_Model.py tests the OOS performance of the in-sample Forward Selection Model
- OLS_Model.py tests the performance of the model with R2 based var selection and OLS coefficient update
- Lasso_Model.py tests the performance of the model with R2 based var selection and Lasso coefficient update
- Fixed_Model*.py tests the performance of the fixed models, there are three parts, covering all the 1-1 pairs from the combined (text and base) var pool
- Stability Checked Model contains all the codes for an alternative fixed model, please refer to the README file there for further details.
#### Codes for Plotting/Interpretation
- 1-1 plot.py plots the time series of the selected variables in the Lasso Update Model.
- best_oneandone.py extracts the winning fixed models for each dependent variable, calculate the MSE ratios against the constant models and save them in .csv
- matrix_plots.py plots the summary matrix for the fixed model. Please refer to the comments in the file for details.


### Notes
- Please change the working directory and the saving directory in each file before running the program.
- Please find the data (v14 for the latest analysis) on the research grid under folder shared/energy_drivers/2020-11-16/data

### Procedures
1. OOS test of the Forward Selection Model
```
chmod u+x OOSfuncs.py

chmod u+x Forward_Model.py

chmod u+x run_Forward_Model.sh

./run_Forward_Model.sh
```
2. OLS Updating Model
```
chmod u+x OOSfuncs.py

chmod u+x OLS_Model.py

chmod u+x run_OLS_Model.sh

./run_OLS_Model.sh
```
3. Lasso Updating Model
```
chmod u+x OOSfuncs.py

chmod u+x Lasso_Model.py

chmod u+x run_Lasso_Model.sh

./run_Lasso_Model.sh
```
4. Fixed Model 
```
chmod u+x Fixed*

chmod u+x Lasso_Model.py

chmod u+x run_Fixed*

sge_run --grid_submit=batch ./run_Fixed*
```
5. Subsequent Plots/Interpretation for the Lasso Model
```
chmod u+x 1-1 plot.py

sge_run --grid_mem=10G --grid_submit=batch './1-1 plot.py'
```
6. Subsequent Plots/Interpretation for the Fixed Model
```
chmod u+x best_oneandone.py

sge_run --grid_mem=10G --grid_submit=batch './best_oneandone.py'

chmod u+x matrix_plots.py

sge_run --grid_mem=10G --grid_submit=batch './matrix_plots.py'
```
