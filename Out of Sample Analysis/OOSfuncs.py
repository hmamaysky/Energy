#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:59:34 2020

@author: Hongyu Wu

This file defines all the functions used in the OOS analysis.
The functions varies from data reader to model selector, all commented
in detail. The main function in each specific code will call them, some 
would be modified a little. 
"""

# %% 0. Importing Pachages

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
import os
import sys
import concurrent.futures
import warnings
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA
warnings.filterwarnings('ignore')

# %% 1. All the functions for the OOS Analysis

### 1.1 Calculate RMSE ###
def RMSE(x):
    ## x is a list recording the difference between the predicted and the real data
    x_s = pd.Series(x).dropna()
    length = len(x_s)
    
    ## return nan if the data is insufficient
    if (length == 0)| (length == 1):
        return np.nan
    
    ## calculating RMSE
    RMSE = np.sqrt(np.sum(np.power(x_s,2))/(length-1))
    return RMSE

### 1.2 Load dataset according to different dependent variable ###
def data_set(d_var):
    ## change directory
    wkdir = '/user/hw2676/files/Energy'
    os.chdir(wkdir)
    
    ## Price vars and Physical vars use different Dataset
    prices_var = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol']

    if d_var in prices_var:
        data = pd.read_stata('data/transformed_data_prices_v12_v4.dta')
        SDFpremium_growing = pd.read_excel('data/SDFgrowing_fut_thurs.xls')
        SDFpremium_rolling = pd.read_excel('data/SDF756rolling_fut_thurs.xls')
    else:
        data = pd.read_stata('data/transformed_data_physical_v12_v4.dta')
        SDFpremium_growing = pd.read_excel('data/SDFgrowing_fut_tues.xls')
        SDFpremium_rolling = pd.read_excel('data/SDF756rolling_fut_tues.xls')
        
    data = pd.merge(data, SDFpremium_growing, on='date', how='left')
    data = pd.merge(data, SDFpremium_rolling, on='date', how='left')
    
    ## Constructed the Sent var here
    data['sent']=data['sCo']+data['sGom']+data['sEnv']+data['sEpg']+data['sBbl']+data['sRpc']+data['sEp']
    
    ## Drop the top 4 rows because no text vars available there
    return data.iloc[4:,:].reset_index(drop=True)

### 1.3 Get the row range of updating and testing data at each window ###
def get_test_row_range(date_col, test_wk, wk=8, update_window=5):
    '''
    Inputs:
        1. date_col: the column with date from the whole dataset
        2. test_wk: the week we are performing algorithm on
        3. wk: 8 indicates LHS is 8wk variable
        4. update_window: years of backward looking for var selection or coefficients update
    
    Outputs:
        1. date_row_update_flag: training rows
        2. date_test_row_flag: testing (forecasting) rows
        3. date_row_pca: rows to calculate first PC using the 5 year window
    '''
    ### Training Window 
    # Here, rows are selected so that the first 8-wk or 4-wk lag does not exceed the updating window
    # Flag 1: Get all the rows with time less or equal to current date
    # Flag 2: Get all rows with lagged time greater than lower bound of the update window
    date_row_flag_1 = (date_col)<=test_wk
    date_row_flag_2 = date_col>(test_wk-pd.Timedelta(str(update_window*365)+'days')+pd.Timedelta(str(7*wk)+ 'days'))
    date_row_update_flag = date_row_flag_1 & date_row_flag_2
    
    ### Testing Window
    # Simply 8 weeks after the current date
    date_test_row_flag = (date_col)==test_wk+pd.Timedelta(str(7*wk)+'days')
    
    ### PCA Window
    # The calculation of the PCA series covers the traning window and 
    # all data within the 8-wk window after the current date to get the PCA values in the test week.
    # Flag 1: Get all rows with lagged time greater than lower bound of the update window
    # Flag 2: Get all the rows with time less or equal to teh testing date
    date_row_flag_1 = date_col>(test_wk-pd.Timedelta(str(update_window*365)+'days')+pd.Timedelta(str(7*wk)+ 'days'))
    date_row_flag_2 = (date_col)<=test_wk+pd.Timedelta(str(7*wk)+'days')
    date_row_pca = date_row_flag_1 &date_row_flag_2
    
    ### Return All the Flags
    return date_row_update_flag, date_test_row_flag, date_row_pca

### 1.4 Get the RHS variable list according to the LHS ###
def ind_var_list(d_var,weeks):
    """
    This function gives the var list on the RHS for a LHS var.
    The list is used to check whether a var is valid to be placed on the RHS,
    as well as pull all the relevant data from the main dataset.
    """
    ### All candidates of RHS vars
    full_list=['DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPIyoy', 'trend', 'VIX', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling',
              'artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']
    ### Price vars will not be included on the RHS for physical dependent vars
    price_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol']
    ### Add d_var if not in RHS list
    ### Here, we use 4wk var for 8wk var training, so no need to add 8wk version in RHS set
    if d_var in price_list:
        full_list.extend(price_list)
    elif d_var not in full_list:
        full_list.insert(0,d_var)
    
    full_list=list(set(full_list))

    return full_list

### 1.5 Dynamically calculate the 1st PCA components in the 5-yr lookback window ###
def PCA_augment(data):
    '''
    This function augments the PCA series for a given lookback window
    '''
    ## Drop fCo in textual var list and freq var list to avoid collinearity
    textual_vars_drop = ['artcount', 'entropy', 'sCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']  
    freq_vars_drop = ['fGom', 'fEnv', 'fEpg', 'fBbl', 'fRpc', 'fEp']
    sent_vars = ['sCo', 'sGom', 'sEnv', 'sEpg', 'sBbl', 'sRpc', 'sEp']
    ## PCA instance with first components on
    pca = PCA(n_components=1)
    data_temp = data.copy()
    ## Augment the data with three PCA series
    data_temp['PCAsent'] = pca.fit_transform(data_temp[sent_vars])
    data_temp['PCAall'] = pca.fit_transform(data_temp[textual_vars_drop])
    data_temp['PCAfreq'] = pca.fit_transform(data_temp[freq_vars_drop])
    ## Return the results
    return data_temp

### 1.6 Select the best vars based on R2 given data
def select_best_n(variables, data_x, y, lags, n=2):
    '''
    Inputs:
        1. variables: a list of baseline or textual variables
        2. data_x: all the RHS data within the training window
        3. y: all the LHS data within the training window
        4. lags: maximum lag in Newey-West robust standard error
    Outputs:
        1. var_1: variable with highest R2
        2. var_2: variable with second highest R2
    '''
    results={}
    vari=[]
    for var in variables:
        x=data_x.loc[:,var]
        
        ### Add a constant column to RHS and regress
        X=sm.add_constant(x, has_constant='add')
        model=lm.OLS(y, X, missing='drop')
        reg=model.fit(cov_type='HAC', cov_kwds={'maxlags':lags})
        
        ### Store the R2 for each variable
        results[var]=reg.rsquared
    ### Get variables with first and second highest R2  
    for i in range(n):
        vari.append(max(results, key=results.get))
        del results[vari[-1]]

    return vari

### 1.7 Select significant vars for base on each 5-yr lookback window
def select_significant(d_var, forecast_start, wk=8, window=5, ns=1):
    '''
    Inputs:
        1. d_var: the LHS variable
        2. forecast_start: date of the current update-forecast window
        3. wk: 8 indicates LHS is 8wk variable
        4. window: years of backward looking for var selection or coefficients update,
                   passed to functions to determine which data are used
    Outputs:
        1. selected: list of vars in the full model
        2. selected_non_textual: list of vars in the baseline model
    '''
    ### Load data and generate lists for baseline and textual variables
    ### Load data and generate lists for baseline and textual variables
    data=data_set(d_var)
    ind_vars=ind_var_list(d_var, weeks=wk)
    if forecast_start<pd.Timestamp('2012-05-11'):
        ind_vars.remove('ovx_cl1')
    if forecast_start<pd.Timestamp('2007-04-06'):
        ind_vars.remove('RPsdf_rolling')
        ind_vars.remove('RPsdf_growing')
    
    textual_vars = ['artcount', 'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']  
    base_vars = [x for x in ind_vars if x not in textual_vars]
    base_vars.remove('sent')

    
    ### Get the training date range
    date_row_range,_,date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    data_y=data[date_row_range]
    
    ### Lag RHS to match LHS data in regression
    # Remember when we choose rows according to update window, we ensure that
    # the lagged observations do not exceed the start line. So for regression,
    # we should move each RHS value 8(4) wk forward to match it with the LHS one.
    lag_vars = ind_vars.copy()
    lag_vars.remove('trend')
    lag_vars.remove('WIPIyoy')
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)
        lag=8

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
        lag=4
        
    data_x_pca=data_x[date_pca_range].loc[:,ind_vars]    
    data_x_pca = PCA_augment(data_x_pca)
    data_x=data_x_pca.iloc[:np.sum(date_row_range),:]
    
    y=data_y[d_var]
    
    textual_vars.extend(['sent', 'PCAfreq','PCAsent', 'PCAall'])
    
    ### Select best four vars based on R-squared
    base_var=select_best_n(base_vars, data_x, y, lags=lag, n=ns)
    text_var=select_best_n(textual_vars, data_x, y, lags=lag, n=ns)
    
    ### Full model and baseline model
    selected = base_var.copy()
    selected.extend(text_var)
    selected_non_textual = base_var
    selected_textual = text_var.copy()
        
    return selected, selected_non_textual, selected_textual