#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:59:34 2020

@author: Hongyu Wu

This file tests the performance of the fixed model. 
A fixed model is a model that the specifications are predetermined during the whole testing horizon.
The 1-1 fixed model takes all the possible combinations from the list of textual and baseline vars.
There are three parts, this is part one, with 1-1 pair indicating a text + a baseline var.
This part kicks off a process for each baseline variables and evaluates the model fixing the baseline var 
as the first fixed model var and setting the second var in turn using all the text variables.
"""
# %% 0. Importing Pachages
import pandas as pd
import sys
import concurrent.futures
import pickle
import os
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
from OOSfuncs import *
__out_dir__ = '/user/hw2676/files/Energy/outputs/wipimom_updated/final_codes_test/fixed_model/'
# %% 1. Defining Functions
### See OOSfuncs.py for all the functions ###

# %% 2. Main Process
def main(var):
    """
    This main function carry out all the steps of the fixed model algorithm:
        1. At each window, look backward for 5 years to update coefficients
        2. Forecast next month or two month (forecast horizon can be controlled in the "get_test_row" function)
        3. Calculate discrepancy between prediction and observation
        4. Calculate and save RMSE for each model and each LHS variables
    Note:
        1. The output is a dataFrames, with 8 dependent variables in the columns and the fixed model specifications in the rows
    """
    ### 0. Set Parameters
    # Read Sys Args
    # weeks=int(sys.argv[1])
    # frequency=int(sys.argv[2])
    # no_varibles=int(sys.argv[3])
    # cv=int(sys.argv[4])
    # base_var=sys.argv[5]
    # text_var=sys.argv[6]

    weeks=8
    frequency=1
    no_varibles=1
    cv=10
    base_var="tnote_10y"
    text_var='artcount'
    # d_var list
    d_vars = [var]
    # Get all coefs for d_var with base_var and text_var
    # These are calculated in the stable model
    os.chdir(__out_dir__)
    coeffs = pickle.load(open('time_varying_model/processed/integrated_var_coefs.p','rb'))[var]
    model = base_var+', '+text_var
    if model in coeffs.keys():
        model_coef = coeffs[model]
    else:
        model_coef = coeffs[text_var+', '+base_var]
        
    # Create lists for result storing
    final_rmse=[]
    # Lookback window (/yr)
    updating_window = 5
    # all the text vars are set to be the second variable in the fixed model in this part
        
    # Result holder: 
    # a dictionary with each dependent variable as key,
    # and a list of rmses for each model as value
    rmse_result={}
    
    ### 1. Main Algorithm
    for d_var in d_vars:
        print('-------------------')
        print(d_var)
        
        # 1.1 Result holder for the prediction differences
        #     Model names are suggested by the var name
        constant_diff_list =[]
        full_diff_list ={}       
        full_diff_list[base_var+', '+text_var]=[]        

        # 1.2 Get the exact date of each window (monthly)
        #     Here, we ensure everytime we train, there are enough observations for backward looking
        #     Note that if a model has ovx_cl1 or sdf variables, the starting time point should be modified accordingly
        time_col = data_set(d_var)['date'] 
        time_lower = time_col[0]
        # Modifying starting date accordingly
        if base_var == 'ovx_cl1':
            time_lower = pd.Timestamp('2007-05-11')
        elif (base_var == 'RPsdf_growing') | (base_var == 'RPsdf_rolling'):
            time_lower = pd.Timestamp('2002-04-08')
        # time_lower = pd.Timestamp('2014-12-12')
        test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
        
        # 1.3 Main algorithm
        # First, at each window, update coefficients of each specifications (model)
        # Second, forecast and get its difference from real observations
        # Finally, record them in a list
        data = data_set(d_var)  # read the data
        ## if the baseline var is a valid RHS var for the dependent variable
        if base_var in ind_var_list(d_var,weeks):
            for week in test_week_list:
                print(week)

                # Update, Forecast and get the prediction error
                full_diff, constant_diff = prediction_one_week(data, d_var, model, model_coef, week)
                # Save the result of the fixed model to corresponding lists
                full_diff_list[model].extend(full_diff)                             
                # Save error of const model to list             
                constant_diff_list.extend(constant_diff)                

            # Store RMSE of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(constant_diff_list),RMSE(full_diff_list[model])]
        ## if the baseline var is not a proper RHS var, just save the constant model results
        else:
            for week in test_week_list:
                print(week)
                
                # Calculate MSE and RMSE
                constant_diff = rolling_diff_Lasso(d_var, [], week, wk=weeks, window=updating_window, cvs=cv)           
                
                # Append results to lists              
                constant_diff_list.extend(constant_diff)                

            # Store RMSE of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(constant_diff_list),RMSE(full_diff_list[model])]

    # 1.4 Create a DataFrame to save results, more convenient for further interpretation
    rmse_df = pd.DataFrame(rmse_result)
    index = ['Constant']
    index.extend([model])
    rmse_df.index=index
    final_rmse.append(rmse_df)

    return final_rmse[0]

# # %% 3. Main Process
# if __name__ == '__main__':
    
#     ### 3.1 System Args    
#     forecasting_week=int(sys.argv[1])
#     update_frequency=int(sys.argv[2])
#     no_variables=int(sys.argv[3])
#     cv_fold=int(sys.argv[4])
#     base_var=sys.argv[5]
#     text_var=sys.argv[6]

#     # forecasting_week=8
#     # update_frequency=1
#     # no_variables=1
#     # cv_fold=10
#     # base_var="WIPImom_8wk"
#     # text_var='sCo'
#     ### 3.2 File Naming Strings
#     if      forecasting_week == 4: file_suffix = '4wk'
#     else:   file_suffix = '8wk'
        
#     if      update_frequency == 1: file_start = 'weekly'
#     else:   file_start = 'monthly'
        
#     if      no_variables  == 1: file = 'oneandone'
#     elif    no_variables == 2: file = 'twoandtwo'
#     else:   file = 'threeandthree'
        
#     if      cv_fold == 5: cv = '5fold'
#     elif    cv_fold == 10: cv = '10fold'
#     else:   cv = '20fold'

#     ### 3.3 Use concurrent calculation for the main algorithm
#     # result holders for each process
#     rmse = pd.DataFrame()
#     # 8 dependent variables
#     d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
#     # Main algorithm and save the results
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         results = executor.map(main, d_vars) 
#         for result in results:
#             rmse = pd.concat([rmse, result], axis=1)
            
#     ### 3.4 Save the results to proper directory
#     rmse.to_excel(__out_dir__+ file+'/'+file_start+'/'+base_var+'/'+base_var+'_Lasso_'+cv+'_'+file_suffix+'_1.xlsx')