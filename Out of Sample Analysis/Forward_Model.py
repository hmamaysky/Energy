#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:36:26 2020

@author: Hongyu Wu

This code use the Lasso coefficient updating method (5-yr lookback window)
to test the out of sample performance of the winning in-sample forward selection models.

The outputs are the MSE ratios of each model and the time series
of the coefficients of the explanatory variables with respect to 
different LHS vars.
"""
# %% 0. Importing Pachages

import pandas as pd
import numpy as np
import concurrent.futures
import statsmodels.api as sm
import sys
import pickle
import warnings
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
from OOSfuncs import *
warnings.filterwarnings('ignore')

# %% 1. Defining Functions for rolling prediction difference (error) calculation
### See OOSfuncs.py for all the functions ###

# %% 2. Main Process
def main(d_var):
    """
    This main function carry out all the steps of the OOS forward model algorithm:
        1. At each window, look backward for 5 years and update coeffs on the forward selection variables.
        2. Forecast and calculate the prediction error
        3. Perform 1.-2. on the Constant Model as well
        4. Calculate and save the RMSEs as well as the coeffs for the forward model
    Outputs:
        1. A DataFrame with RMSEs of the constant and the forward selection model
        2. A dictionary with the time series of the coeffs of each variable in the forward selection model
    """
    ### 0. Sys args for further usage
    weeks=int(sys.argv[1])
    frequency=int(sys.argv[2])
    no_varibles=int(sys.argv[3])
    cv=int(sys.argv[4])
    
    ### 1. The vars of the forward selection models from the in-sample analysis
    forward_selected_models={'FutRet':['FutRet', 'DInv', 'DFX', 'basis', 'WIPIyoy', 'VIX', 'PCAsent', 'DOilVol'],
                             'xomRet':['xomRet', 'sp500Ret', 'WIPIyoy', 'tnote_10y', 'DInv', 'VIX', 'DFX', 'sRpc'],
                             'bpRet':['bpRet', 'sp500Ret', 'sEp', 'DFX', 'DSpot', 'fEp', 'sGom', 'tnote_10y'],
                             'rdsaRet':['rdsaRet', 'fBbl', 'sEnv', 'fGom', 'sEp', 'DInv', 'sCo', 'VIX'],
                             'DSpot':['DSpot', 'fRpc', 'fBbl', 'basis', 'sp500Ret', 'sEnv', 'sGom', 'entropy'],
                             'DOilVol':['DOilVol', 'OilVol', 'DSpot', 'VIX', 'entropy', 'fGom', 'fCo', 'PCAsent'],
                             'DInv':['DInv', 'DProd', 'artcount', 'fRpc', 'WIPIyoy', 'entropy', 'sEp', 'vix_spx'],
                             'DProd':['DProd', 'sEp', 'DInv', 'sBbl', 'fRpc', 'DOilVol', 'WIPIyoy', 'VIX']}

    ### 2. Data Preparation
    # 2.1 Full Dependent Variable List for later loop
    d_vars = [d_var]
    # 2.2 Create lists for result storing
    updating_window = 5

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print("Update Window: ", updating_window)
    
    # 2.3 Result holder: 
    # a dictionary with each dependent variable as key,
    # and a list of rmses for each model as value
    rmse_result={}
    var_coef_result={}
    
    for d_var in d_vars:
        print('-------------------', d_var, sep='\n')
        ### 2.4 Read in data and Some result holders
        # 2.4.1 Read in the whole dataset
        data = data_set(d_var)
        # 2.4.2 Get the specification of the forward model (top 2 or top 7)
        forward_vars = forward_selected_models[d_var][:no_varibles+1]
        # 2.4.3 lists to store the diff between the prediction and the real data
        constant_diff_list = []
        forward_diff_list = []
        # 2.4.4 dict to store the coeff of each var
        var_coeff_ = dict()     
        # 2.4.5 Get the exact date of each window (monthly)
        # Here, we ensure everytime we train, there are enough observations for backward looking
        time_col = data_set(d_var)['date'] 
        test_week_list = [time for time in time_col if time>=time_col[0]+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]

        ### 3. Main algorithm
        for week in test_week_list:
            print(week)
            
            # 3.1 Calculate prediction error for the constant model and the forward model
            #     Store the coeff for the forward model as well
            constant_diff,_ = rolling_diff_forward(data, d_var, [], week, wk=weeks, window=updating_window, cvs=cv)
            forward_diff, var_coeff = rolling_diff_forward(data, d_var, forward_vars, week, wk=weeks, window=updating_window, cvs=cv)
            
            
            # 3.2 Append diffs to lists              
            constant_diff_list.extend(constant_diff)                
            forward_diff_list.extend(forward_diff)
            
            # 3.3 Save the var coeffs
            var_coeff_[week]=var_coeff
            
        # 3.4 Store RMSE of different models for each LHS variable  
        rmse_result[d_var]=[RMSE(constant_diff_list), RMSE(forward_diff_list)]
        var_coef_result[d_var] = pd.DataFrame(var_coeff_)
        var_coef_result[d_var].index = ['Constant']+forward_vars


    ### 4. Process and Return the results
    rmse_df = pd.DataFrame(rmse_result)
    rmse_df.index=['Constant', 'Forward Model']

    return rmse_df, var_coef_result

# %% 3. Main Process
if __name__ == '__main__':
    ### Sys args for file naming
    forecasting_week=int(sys.argv[1])
    no_variables=int(sys.argv[3])

    ### File naming 
    if forecasting_week == 4:
        file_suffix = '4wk'
    else:
        file_suffix = '8wk'
    
    ### Use concurrent calculation for the main algorithm
    # result holders for each process
    rmse_df=pd.DataFrame()
    var_coef_dict = dict()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # 8 dependent vars
        d_var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
        # Main algorithm and save the results
        results = executor.map(main, d_var_list)
        try:
            for result in results:
                rmse_df = pd.concat([rmse_df, result[0]], axis=1)
                var_coef_dict.update(result[1])
        except:
            pass
    # Save the outputs to certain directory
    rmse_df.to_excel('/user/hw2676/files/Energy/outputs/model_selection/parsimonious/forward_models/Lasso_10fold_'
                    +file_suffix+str(no_variables)+'vars_forward.xlsx')
    pickle.dump(var_coef_dict, open('/user/hw2676/files/Energy/outputs/model_selection/parsimonious/forward_models/var_coef_'
                                +str(no_variables)+'vars_forward.p','wb'))
