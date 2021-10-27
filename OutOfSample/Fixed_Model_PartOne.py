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
__out_dir__ = '/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/'
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
    weeks=int(sys.argv[1])
    frequency=int(sys.argv[2])
    no_varibles=int(sys.argv[3])
    cv=int(sys.argv[4])
    base_var=sys.argv[5]

    # weeks=8
    # frequency=1
    # no_varibles=1
    # cv=10
    # base_var="ovx_diff"
    # d_var list
    d_vars = [var]

    # Create lists for result storing
    final_rmse=[]
    # Lookback window (/yr)
    updating_window = 5
    # all the text vars are set to be the second variable in the fixed model in this part
    textual_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAfreq','PCAsent', 'PCAall']
    # textual_vars = ['artcount']
        
    # Result holder: 
    # a dictionary with each dependent variable as key,
    # and a list of rmses for each model as value
    rmse_result={}
    # Get the dictionary saving all the pre-generated coefficients
    os.chdir(__out_dir__)
    coeffs = pickle.load(open('time_varying_model/processed/integrated_var_coefs.p','rb'))[var]
    
    ### 1. Main Algorithm
    for d_var in d_vars:
        print('-------------------')
        print(d_var)
        
        # 1.1 Result holder for the prediction differences
        #     Model names are suggested by the var name
        constant_diff_list =[]
        constant_RMSE_list = []
        full_RMSE_list ={} 
        full_diff_list ={}
        model_coef = {}
        for t_var in textual_vars:
            full_diff_list[base_var+', '+t_var]=[] 
            full_RMSE_list[base_var+', '+t_var]=[]  
            model = base_var+', '+t_var
            if model in coeffs.keys():
                model_coef[model] = coeffs[model]
            else:
                model_coef[model] = coeffs[t_var+', '+base_var]
        
        # 1.2 Get the exact date of each window (monthly)
        #     Here, we ensure everytime we train, there are enough observations for backward looking
        #     Note that if a model has ovx_cl1 or sdf variables, the starting time point should be modified accordingly
        time_col = data_set(d_var)['date'] 
        time_lower = time_col[0]
        # time_lower = pd.Timestamp('2015-02-01')
        # Modifying starting date accordingly
        if base_var == 'ovx_diff':
            time_lower = pd.Timestamp('2007-05-11')
        elif (base_var == 'RPsdf_growing') | (base_var == 'RPsdf_rolling'):
            time_lower = pd.Timestamp('2002-04-08')
        # time_lower = pd.Timestamp('2015-01-02')
        test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
        
        # 1.3 Main algorithm
        # First, at each window, update coefficients of each specifications (model)
        # Second, forecast and get its difference from real observations
        # Finally, record them in a list
        data=data_set(d_var)
        ## if the baseline var is a valid RHS var for the dependent variable
        if base_var in ind_var_list(d_var,weeks):
            for week in test_week_list:
                print(week)
                # Update, Forecast and get the prediction error
                constant_diff = [rolling_diff(data, d_var, [], week)]
                for t_var in textual_vars:
                    print(t_var)
                    model = base_var+', '+t_var
                    full_diff, _ = prediction_one_week(data, d_var, model, model_coef[model], week)
                    # Save the result of the fixed model to corresponding lists
                    full_diff_list[model].extend(full_diff)
                    full_RMSE_list[model].append((RMSE(full_diff_list[model])))
                
                # Save error of const model to list             
                constant_diff_list.extend(constant_diff)                
                constant_RMSE_list.append(RMSE(constant_diff_list))
            
            # Store RMSE of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(constant_diff_list)]
            rmse_result[d_var].extend([RMSE(full_diff_list[base_var+', '+x]) for x in textual_vars])
        ## if the baseline var is not a proper RHS var, just save the constant model results
        else:
            for week in test_week_list:
                print(week)
                
                # Calculate MSE and RMSE
                constant_diff = [rolling_diff(data, d_var, [], week)]         
                
                # Append results to lists              
                constant_diff_list.extend(constant_diff)                

            # Store RMSE of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(constant_diff_list)]
            rmse_result[d_var].extend([RMSE(full_diff_list[base_var+', '+x]) for x in textual_vars])

    # 1.4 Create a DataFrame to save results, more convenient for further interpretation
    rmse_df = pd.DataFrame(rmse_result)
    index = ['Constant']
    index.extend([base_var+', '+x for x in textual_vars])
    rmse_df.index=index
    final_rmse.append(rmse_df)
    
    # 1.5 Save MSE results
    rmse_dict={d_var+'_fullrmse_1':full_RMSE_list,d_var+'_consrmse':constant_RMSE_list}

    return final_rmse[0], rmse_dict

# %% 3. Main Process
if __name__ == '__main__':
    
    ### 3.1 System Args    
    forecasting_week=int(sys.argv[1])
    update_frequency=int(sys.argv[2])
    no_variables=int(sys.argv[3])
    cv_fold=int(sys.argv[4])
    base_var=sys.argv[5]

    # forecasting_week=8
    # update_frequency=1
    # no_variables=1
    # cv_fold=10
    # base_var="HedgPres"
    ### 3.2 File Naming Strings
    if      forecasting_week == 4: file_suffix = '4wk'
    else:   file_suffix = '8wk'
        
    if      update_frequency == 1: file_start = 'weekly'
    else:   file_start = 'monthly'
        
    if      no_variables  == 1: file = 'oneandone'
    elif    no_variables == 2: file = 'twoandtwo'
    else:   file = 'threeandthree'
        
    if      cv_fold == 5: cv = '5fold'
    elif    cv_fold == 10: cv = '10fold'
    else:   cv = '20fold'

    ### 3.3 Use concurrent calculation for the main algorithm
    # result holders for each process
    rmse = pd.DataFrame()
    rmse_dicts={}
    # 8 dependent variables
    d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    # Main algorithm and save the results
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(main, d_vars) 
        for result in results:
            rmse = pd.concat([rmse, result[0]], axis=1)
            rmse_dicts.update(result[1])
            
    ### 3.4 Save the results to proper directory
    rmse.to_excel('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/'
                    +file+'/'+file_start+'/raw/'+base_var+'_Lasso_'+cv+'_'+file_suffix+'_1.xlsx')
    pickle.dump(rmse_dicts,open('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/'
                    +file+'/'+file_start+'/rmse/'+base_var+'_Lasso_'+cv+'_'+file_suffix+'_1.p','wb'))
    