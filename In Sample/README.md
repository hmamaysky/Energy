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
5. Subsequent Processing for the Fixed Model
```
chmod u+x best_oneandone.py

sge_run --grid_mem=10G --grid_submit=batch './best_oneandone.py'

chmod u+x best_ratio.py

sge_run --grid_mem=10G --grid_submit=batch './best_ratio.py'

chmod u+x matrix_plots.py

sge_run --grid_mem=10G --grid_submit=batch './matrix_plots.py'
```

