# Energy
## Out of Sample Analysis

Contains codes for OOS analysis for the Energy Project

Notes:
- OOSfuncs.py includes all the functions needed for the OOS analysis. 
- Forward_Model.py tests the OOS performance of the in-sample Forward Selection Model
- OLS_Model.py tests the performance of the model with R2 based var selection and OLS coefficient update
- Lasso_Model.py tests the performance of the model with R2 based var selection and Lasso coefficient update

Procedures:
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
