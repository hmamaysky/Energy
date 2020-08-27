# Stability Checked Model
## A Variation from the Fixed Model

This fixed model deals with the time variation of the coefficients of the selected models.

A model becomes a prediction candidate on a forecastting date only if it has stable and nonzero coefficients over
the past few years. (being stable means no switch of sign)

To perform the analysis, run the files in following order:

1. shuffle_model_groups.py
    
    This file separate all possible combinations of the 39 (20 text + 19 baseline) RHS variables into 30 groups.
    
    The separation will speed up the next step by leveraging the parallelization techniques.
    
    ```
    chmod u+x OOSfuncs.py
    chmod u+x shuffle_model_groups.py
    sge_run --grid_mem=32G --grid_submit=batch "./shuffle_model_groups.py"
    ```
    
2. run_generate_model_coefs.sh
    
    This file kicks off 30 generate_model_coefs_v1.0.py in parallel, which calculate all the Lasso 5yr lookback coefficients.
    
    The results are saved in 30 separated .p files for further usage.
    
    ```
    chmod u+x generate_model_coefs.py
    chmod u+x run_generate_model_coefs.sh
    ./run_generate_model_coefs.sh
    ```
    
3. run_select_models.sh
    
    This file runs select_models.py, which aggregates all the 30 .p files above and selects the ones that
    qualify as a stable or nonzero model.
    
    The outputs are dictionaries with LHS vars as keys, with some dictionaries as values. 
    
    These dictionaries have weektime as key and another layer of dict storing all the selected models and their coefs as values.

    ```
    chmod u+x select_models.py
    chmod u+x run_select_models.sh
    ./run_select_models.sh
    ```
    
4. run_Time_varying_model.sh
    
    This file calls Time_varying_model.py, which performs prediction weekly using the selected models
    and their coefficients. Constant models of the same period are calculated as well.
    
    Two models are tested here: 1. the plain time varying model, 2. time varying model filled with the constant prediction 
    if there isn't a valid prediction candidate. 
    
    The reported results are presented in an excel file with the RMSE ratios of the constant and time varying model.

    ```
    chmod u+x Time_varying_model.py
    chmod u+x run_Time_varying_model.sh
    ./run_Time_varying_model.sh
    ```
