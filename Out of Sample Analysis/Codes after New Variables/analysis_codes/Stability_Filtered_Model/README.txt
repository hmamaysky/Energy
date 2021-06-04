This fixed model deals with the time variation of selected models.
A model become a prediction candidate only if it has stable and nonzero coefficients over
the past few years. (being stable means no switch of sign)

To perform the analysis, run the files in following order:
1. shuffle_model_groups.py
    This file separate all possible combinations of the 47 (20 text + 19 baseline (+ 8 new)) RHS variables into 30 groups.
    The separation will speed up the next step by leveraging the parallelization techniques.
    
2. run_generate_model_coefs_v1.0.sh
    This file kicks off 30 generate_model_coefs_v1.0.py in parallel, which calculate all the Lasso 5yr lookback coefficients.
    The results are saved in 30 separated .p files for further usage.
    
3. run_select_models.sh
    This file runs select_models.py, which aggregates all the 30 .p files above and selects the ones that
    qualify as a stable or nonzero model.
    The outputs are dictionaries with LHS vars as keys, with some dictionaries as values. 
    These dictionaries have weektime as key and another layer of dict storing all the selected models and their coefs as values.
    
4. run_Time_varying_model.sh
    This file calls Time_varying_model.py, which performs prediction weekly using the selected models
    and their coefficients. Constant models of the same period are calculated as well.
    Two models are tested here: 1. the plain time varying model, 2. time varying model filled with the constant prediction 
    if there isn't a valid prediction candidate. 
    The reported results are presented in an excel file with the RMSE ratios of the constant and time varying model.
