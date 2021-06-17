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
as the first fixed model var and setting the second var in turn using all the baseline variables.
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
    
    # d_var list
    d_vars = [var]

    # Create lists for result storing
    final_rmse=[]
    updating_window = 5
    # all the var combinations are based on the baseline var list in this part
    textual_vars = ["FutRet", "xomRet", "bpRet", "rdsaRet", "DOilVol",
                    "OilVol", "DInv", "DProd", "DSpot", "tnote_10y",
                    "DFX", "sp500Ret", "basis", "WIPImom_{}wk".format(weeks), "trend", 'bmratio', 'mom_fut1', 'basismom', 'dxy_betas', 'cpiyr_betas', 'hp', 'liquidity', 'oi',
                    "RPsdf_growing", "RPsdf_rolling", "vix_spx", "ovx_cl1"]
    # remove the first var from the list for the second var to avoid invalid model
    try:
        textual_vars.remove(base_var)
    except:
        pass
    
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
        
        # 1.1 Make sure all the vars in the list are valid for the given LHS dependent variable
        textual_vars = [x for x in textual_vars if x in ind_var_list(d_var,weeks)]
        
        # 1.2 Result holder for the prediction differences
        #     Model names are suggested by the var name
        full_diff_list ={}
        model_coef = {}
        model_dict = {} # coef has order, this is crucial in prediction process (multiply the correct coef)
        for t_var in textual_vars:
            full_diff_list[base_var+', '+t_var]=[]   
            model = base_var+', '+t_var
            if model in coeffs.keys():
                model_coef[model] = coeffs[model]
                model_dict[model] = model
            else:
                model_coef[model] = coeffs[t_var+', '+base_var]  
                model_dict[model] = t_var+', '+base_var
        # 1.3 Get the exact date of each window (monthly)
        #     Here, we ensure everytime we train, there are enough observations for backward looking
        #     Note that if a model has ovx_cl1 or sdf variables, the starting time point should be modified accordingly
        time_col = data_set(d_var)['date']         
        
        # 1.4 Main algorithm
        # First, at each window, update coefficients of each specifications (model)
        # Second, forecast and get its difference from real observations
        # Finally, record them in a list
        data=data_set(d_var)
        for t_var in textual_vars:
            # Modifying starting date accordingly
            time_lower = time_col[0]
            if (base_var == 'RPsdf_growing') | (base_var == 'RPsdf_rolling') | (t_var == 'RPsdf_growing') | (t_var == 'RPsdf_rolling'):
                time_lower = pd.Timestamp('2002-04-08')
            elif (base_var =='ovx_cl1') |(t_var == 'ovx_cl1'):
                time_lower = pd.Timestamp('2007-05-11')
            # Get time horizon first, because the second element in a pair may have different time horizon than others
            test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
            
            # Update, Forecast and get the prediction error
            if base_var in ind_var_list(d_var,weeks):
                for week in test_week_list:
                    print(week)
                    model = base_var+', '+t_var
                    full_diff, _ = prediction_one_week(data, d_var, model_dict[model], model_coef[model], week)
                    full_diff_list[base_var+', '+t_var].extend(full_diff)
        # 1.5 Save the results
        rmse_result[d_var]=[]
        rmse_result[d_var].extend([RMSE(full_diff_list[base_var+', '+x]) for x in textual_vars])
        
    # 2. Create a DataFrame to save results, more convenient for further interpretation
    rmse_df = pd.DataFrame(rmse_result)
    index = []
    index.extend([base_var+', '+x for x in textual_vars])
    rmse_df.index=index
    final_rmse.append(rmse_df)

    return final_rmse[0]

# %% 3. Main Process
if __name__ == '__main__':
    
    ### 3.1 System Args    
    forecasting_week=int(sys.argv[1])
    update_frequency=int(sys.argv[2])
    no_variables=int(sys.argv[3])
    cv_fold=int(sys.argv[4])
    base_var=sys.argv[5]

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
    # 8 dependent variables
    d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    # Main algorithm and save the results
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(main, d_vars) 
        for result in results:
            rmse = pd.concat([rmse, result], axis=1)
            
    ### 3.4 Save the results to proper directory
    rmse.to_excel('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/'
                    +file+'/'+file_start+'/'+base_var+'_Lasso_'+cv+'_'+file_suffix+'_3.xlsx')