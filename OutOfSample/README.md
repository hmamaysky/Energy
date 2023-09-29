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


### Notes before Running Analysis
- Please change the working directory and the saving directory in each file before running the program.
- Please find the data (v14 for the latest analysis) on the research grid under folder shared/energy_drivers/2020-11-16/data

### Pre-steps
1. Obtain rolling topic models using Louvain and Genetic Algorithm, and store the membership

2. Generate clustering_C.csv

3. Run topic allocation for each month based on backward-looking topic models (less than 3min)
```
./topic_allocation.py --inputWordsPath='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_clustering_C' --outputPath='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_topic_allocation' --local_topic_model=True
```

4. Combine info (less than 3min)
```
./info.py --outputPath='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info' --local_topic_model=True
```
   

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

Here, we run all the steps in Stability Checked model first to save time on data processing. Please refer to the folder for details. After that, run following codes.
```
chmod u+x Fixed*

chmod u+x Lasso_Model.py

chmod u+x run_Fixed*

sge_run --grid_submit=batch ./run_Fixed*
```

