#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:41:54 2020

@author: Hongyu Wu

Description:
------------
    This file takes in the calculated var coefs (from generate_model_coef.py) and aggregate them into one file.
    Then selects the stable models and saves them in a dictionary for each week.
    By stable, it means the coefs do not hit zero or change sign in the past several years 
    (specified by the researcher, see var 'lookback' in the main process). 
Note:
-----    
    Please modify the directories in the "read_and_integrate" function, the "main" function
    and the main process to generate desired results.
"""

# %% 0. Importing Packages
import pandas as pd
import os
import pickle
import concurrent.futures
import functools
import sys
from OOSfuncs import *

# %% 1. Defining Functions
def read_and_integrate(first_time=0):
    '''
    If first_time=0:
        Read all the 30 pickle files with 8 LHS var as keys and dictionaries as values.
        The values are dictionaries as well, with model names as keys and coef series as values.
        Then integrate them, save to local, and return.
    else:
        read the preciously integrated and saved file and return it
    '''
    if first_time:
        files=os.listdir()
        whole_dict=pickle.load(open(files[0],'rb'))
        for i in range(1,len(files)):
            partial_dict=pickle.load(open(files[i],'rb'))
            for key in whole_dict.keys():
                whole_dict[key].update(partial_dict[key])
        pickle.dump(whole_dict, open('../processed/integrated_var_coefs.p','wb'))
    else:
        whole_dict=pickle.load(open('../processed/integrated_var_coefs.p','rb'))
    return whole_dict

def select_stable_one_week(coef_df, timestamp, lookback=3):
    '''
    Check if nonzero or stable for a single model at a weekly checkpoint.
    Return the last coefs if past, 0 if not.
    '''
    # Time stamps in the lookback window
    times = [timestamp-pd.Timedelta(str(7*(lookback*52-i))+'days') for i in range(lookback*52+1)]
    # All coefs within time window
    try:
        coefs = coef_df.loc[:,times]
    except:
        return pd.Series()
    # The first and second coef series (Constant is the 0th, ignore it when testing, but need it if qualified) 
    coef_a = coefs.iloc[1,:]
    coef_b = coefs.iloc[2,:]
    # First criterion: Nonzero throughout the period
    a_nonzero = all(x!=0 for x in coef_a)
    b_nonzero = all(x!=0 for x in coef_b)
    # Second criterion: No switching of sign 
    a_stable = max(all(x>0 for x in coef_a),all(x<0 for x in coef_a))
    b_stable = max(all(x>0 for x in coef_b),all(x<0 for x in coef_b))
    # If all passed, the min value should be 1, then return the coefs and contant at the checkpoint
    # These coefs will be used for further prediction
    # If not qualified, return 0 so the model won't enter the final dictionary
    if min(a_nonzero, b_nonzero, a_stable, b_stable):
        return coefs.iloc[:,-1]
    else:
        return pd.Series()
    

def select_stable_whole_sample(d_var, var_coefs_dict, lookback=3):
    '''
    For a given LHS var, and for each week, select the stable models within the lookback window.
    Return the dictionary with timestamps as keys and selected models with lastest coefficients as values.
    Inputs:
        1. d_var: One of the eight LHS vars
        2. vars_coefs_dict: the integrated var coef dict, get from the read_and_intergrate function
        3. lookback: look back window for stability checking
    Output:
        1. a dictionary with timestamp as keys and dicts as values where the dicts has selected
           models as keys and the latest coefs as values.
    '''
    # The result dict with qualified models and their coefs for each week
    result_df=dict()
    result_df[d_var]=dict()
    # Read the coef dictionary for a particular LHS var
    var_coefs = var_coefs_dict[d_var]
    # Separate it because ovx and sdf related models have shorter time span
    var_coefs_ovx = {key:var_coefs[key] for key in var_coefs.keys() if 'ovx' in key}
    var_coefs_sdf = {key:var_coefs[key] for key in var_coefs.keys() if ('sdf' in key) and (key not in var_coefs_ovx.keys())}
    var_coefs = {key:var_coefs[key] for key in var_coefs.keys() if (key not in var_coefs_ovx.keys()) and (key not in var_coefs_sdf.keys())}
    # Get time horizon
    time_col = data_set(d_var)['date']
    time_lower = time_col[0]
    week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*(5+lookback)*52)+'days')][::1]
    # Perform selection
    for week_time in week_list:
        print(week_time)
        # Handle ovx and sdf related models by one time update at the startig point
        if week_time==pd.Timestamp('2007-05-11')+pd.Timedelta(str(7*(5+lookback)*52)+'days'):
            var_coefs.update(var_coefs_ovx)
        if week_time==pd.Timestamp('2002-04-08')+pd.Timedelta(str(7*(5+lookback)*52)+'days'):
            var_coefs.update(var_coefs_sdf)
        # create a dict for a certain week
        result_df[d_var][week_time]=dict()
        # core process for selection
        ## recall that the keys of the dict are model names, and values are coef series.
        for model,coefs in var_coefs.items():
            # select result for that week
            select_result = select_stable_one_week(coefs, week_time, lookback=lookback)
            # if selection is not empty, add it to the weekly dictionary
            if any(select_result): result_df[d_var][week_time][model]=select_result
    
    return result_df
    
    
# %% 2. Main Function    
def main(d_var, lookback=3):
    # Set working directory
    wkdir = '/user/hw2676/files/Energy/outputs/model_selection/fixed_model/time_varying_model/raw'
    os.chdir(wkdir)
    # Read and integrate all the models, 
    # The integrated results will be saved if first_time is set to be 1
    # Later change first_time to 0 to avoid redundant aggregation process
    var_coefs_dict = read_and_integrate(first_time=1)
    # Generate all the stable models for the whole sample period
    selected_models = select_stable_whole_sample(d_var, var_coefs_dict, lookback)
    
    return selected_models

# %% 3. Main Process
if __name__=='__main__':
    # Get lookback window from command line
    lookback = int(sys.argv[1])
    # Parallelize the whole process to speed up the calculation
    dvar_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    selected_models_dict = dict()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(functools.partial(main, lookback=lookback), dvar_list)
        for result in results:
            selected_models_dict.update(result)
    # Save the results to proper directory
    pickle.dump(selected_models_dict, 
                open('/user/hw2676/files/Energy/outputs/model_selection/fixed_model/'+
                     'time_varying_model/processed/'+str(lookback)+'yrback_selected_models_with_coefs.p',
                     'wb'))
        
        
        