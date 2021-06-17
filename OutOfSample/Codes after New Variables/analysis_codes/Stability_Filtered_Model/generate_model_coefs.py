#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 20:45:04 2020

@author: Hongyu Wu

Description: 
------------
    This file kicks off 30 generate_model_coefs_v1.0.py in parallel, which calculate all the Lasso 5yr lookback coefficients.
    The results are saved in 30 separated .p files for further usage.
Note:
-----
    Please change the directories in the "main" function and the main process accordingly before running the program.
"""


# %% 0. Importing Pachages
import pandas as pd
import concurrent.futures
import sys
import pickle
from collections import OrderedDict
from OOSfuncs import *

# %% 1. Defining Functions
### See OOSfuncs.py for all the functions ###

# %% 2. Main Function
def main(d_var):
    """
    This main function carry out all the steps of coefficient calculation:
        1. At each window, look backward for 5 years and select RHS variables.
        2. Update coefficients of the selected ones using that 5 year window again.
        3. Forecast next month or two month (this can be controlled in the "get_test_row" function)
        4. Calculate coefficients for each model on each prediction date
        5. Save the results in a dictionary
    Note:
        1. A dictionary, with d_var as key, a dictionary as value
           The second dictionary has model specifications ("var1, var2") as keys, 
           and time series of the coefficients as values (saved in a pd.DataFrame).
    """
    ### 0. Set Parameters
    # Read Sys Args
    weeks=int(sys.argv[1])
    frequency=int(sys.argv[2])
    cv=int(sys.argv[3])
    n_round=sys.argv[4]

    # weeks=8
    # frequency=1
    # cv=10
    # n_round='30' 
    # Read in all the models for round n
    model_list = pickle.load(open('shuffled_model_list.p','rb'))[n_round]
    
    # d_var list
    d_vars = [d_var]

    # Lookback window (/yr)
    updating_windows = [5]
    
    # Dictionary to save coefs series for a specific model
    var_coef_dict = OrderedDict()
    
    # Set the key and value (another dict, date as key, coefs as value)
    for model in model_list:
        var_coef_dict[model[0]+', '+model[1]]=OrderedDict()
        
    # Check if the model has invalid RHS var for the d_var
    valid_ind_vars = ind_var_list(d_var,8) + ['PCAsent', 'PCAfreq', 'PCAall']
    model_list = [model for model in model_list if (model[0] in valid_ind_vars) and (model[1] in valid_ind_vars)]
    print('*************************')
    print('*************************')

    ### 1. Main Algorithm
    for updating_window in updating_windows:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print("Update Window: ", updating_window)
        
        for d_var in d_vars:
            print('-------------------')
            print(d_var)
            
            # 1.1 Get the coeff for each model using Lasso regression on 5-yr lookback windows
            for model in model_list: 
                ### Read in dataset
                data = data_set(d_var)
                ### Set time lower bound and get all the forecasting weeks
                # if model has ovx_cl1 or SDF related vars, the starting time is postponed
                # due to data issue
                time_col = data_set(d_var)['date']
                time_lower = time_col[0]
                if 'ovx_cl1' in str(model):
                    time_lower = pd.Timestamp('2007-05-11')
                elif ('RPsdf_growing' in str(model)) | ('RPsdf_rolling' in str(model)):
                    time_lower = pd.Timestamp('2002-04-08')
                test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
                # Get the coefficients for each week and save them
                for week in test_week_list:
                    print(week)
                    _, var_coeff = rolling_diff_stability_coef(data, d_var, model, week, wk=weeks, window=updating_window, cvs=cv)
                    var_coef_dict[model[0]+', '+model[1]][week]=var_coeff
        # 1.2 Save the results
        ### Each key corresponds to a DF with time points as columns and Intercept, var1, var2, R2 as row names
        for model in model_list:
            var_coef_dict[model[0]+', '+model[1]] = pd.DataFrame(var_coef_dict[model[0]+', '+model[1]],
                         index=['Constant',model[0],model[1],'R2'])
    ### 2. Return the results
    # Save in a dict for later aggregation
    d_var_dict = dict()
    d_var_dict[d_var] = var_coef_dict
    return d_var_dict

# %% 3. Main Process
if __name__ == '__main__':

    ### 3.1 System Args
    forecasting_week=int(sys.argv[1])
    update_frequency=int(sys.argv[2])
    cv_fold=int(sys.argv[3])
    n_rounds = sys.argv[4]

    # forecasting_week=8
    # update_frequency=1
    # cv_fold=10
    # n_rounds = '30'
    ### 3.2 File Naming Strings
    if forecasting_week == 4: file_suffix = '4wk'
    else: file_suffix = '8wk'
                
    if cv_fold == 5: cv = '5fold'
    elif cv_fold == 10: cv = '10fold'
    else: cv = '20fold'
    
    ### 3.3 Use concurrent calculation for the main algorithm
    # result holders for each process
    var_coef_dict = dict()
    # 8 dependent variables
    d_var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    # d_var_list = ['DInv', 'DProd']
    with concurrent.futures.ProcessPoolExecutor() as executor:     
        results = executor.map(main, d_var_list)
        # try:
        for result in results:
            var_coef_dict.update(result)
        # except:
        #     pass
            
    ### 3.4 Save the results to proper directory     
    pickle.dump(var_coef_dict, open('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/raw/var_coefs_round_'
                                +n_rounds+'.p','wb'))