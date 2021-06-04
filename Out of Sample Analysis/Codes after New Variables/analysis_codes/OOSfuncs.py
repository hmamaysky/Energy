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
        data = pd.read_stata('data/transformed_data_prices_v17.dta').rename(columns={'date_Fri':'date'})
        SDFpremium_growing = pd.read_excel('data/SDFgrowing_fut_thurs.xls')
        SDFpremium_rolling = pd.read_excel('data/SDF756rolling_fut_thurs.xls')
    else:
        data = pd.read_stata('data/transformed_data_physical_v17.dta').rename(columns={'date_Tue':'date'})
        data['date'] = data['date'].apply(lambda x:x+pd.Timedelta('3 days'))
        SDFpremium_growing = pd.read_excel('data/SDFgrowing_fut_tues.xls')
        SDFpremium_rolling = pd.read_excel('data/SDF756rolling_fut_tues.xls')
    
    # Rename columns, expect the date variables and the trend
    all_columns = list(data.columns.values)
    all_columns.remove('trend')
    all_columns = [x for x in all_columns if 'date' not in x]
    
    data = data.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(all_columns)})
    data = pd.merge(data, SDFpremium_growing, on='date', how='left')
    data = pd.merge(data, SDFpremium_rolling, on='date', how='left')
    
    ## Constructed the Sent var here
    data['sent']=data['sCo']+data['sGom']+data['sEnv']+data['sEpg']+data['sBbl']+data['sRpc']+data['sEp']
    
    ## Drop the top 4 rows because no text vars available there
    return data.iloc[3:,:].reset_index(drop=True)

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
    full_list=['DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot', 'FutRet', 'xomRet', 'bpRet', 'rdsaRet',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPImom_{}wk'.format(weeks), 'trend', 'VIX', 'vix_diff', 'ovx_diff', 'RPsdf_growing', 'RPsdf_rolling',
              'BEME', 'Mom', 'BasMom', 'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt',
              'artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']
    ### Price vars will not be included on the RHS for physical dependent vars
    price_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol']
    # price_list = []
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
        ind_vars.remove('ovx_diff')
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
    lag_vars.remove('WIPImom_{}wk'.format(wk))
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

# %% 2. Main algorithms for each week
### 2.1 Custom Rolling Diff function for OLS Updating method
def rolling_diff_OLS(d_var, ind_vars, forecast_start, wk=8, window=5):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_start: current week which we base on to forecast 
        4. wk: 8 or 4 according to which vars we are interested in
        5. zero: True if 0 Specification model
        6. window: The backward looking length for coefficient updating
    Outputs:
        1. diff: difference between forecast and real observation
    '''

    ### We need to update and then forecast
    data=data_set(d_var)
    
    ### Get update and forecast window
    date_update_range,date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    
    ### Shift x to match y and set Newey-West max lag
    lag_vars = ind_var_list(d_var, weeks=wk)
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_{}wk'.format(wk))
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)
        lag=8

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
        lag=4
    
    ### train suffix is the update window, test suffix is the forecast window
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### Update PCA values        
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)
    
    data_xtrain= data_x_pca.iloc[:np.sum(date_update_range),:]
    data_xtest = data_x_pca.iloc[[-1],:]   
    
    ### Run regression and update, here, sm.add_constant adds a constant to RHS
    x_train=data_xtrain.loc[:,ind_vars]
    X_train=sm.add_constant(x_train)
    if wk==8:
        y_train=data_ytrain[d_var+'_t8']
    else:
        y_train=data_ytrain[d_var+'_t4']
    model=lm.OLS(y_train, X_train, missing='drop')
    reg=model.fit(cov_type='HAC', cov_kwds={'maxlags':lag})
    x_test=data_xtest.loc[:,ind_vars]
    X_test=sm.add_constant(x_test, has_constant='add')
    
    ### Predict and record the difference
    if wk==8:
        y_test=data_ytest[d_var+'_t8']
    else:
        y_test=data_ytest[d_var+'_t4']
    diff = reg.predict(X_test) - y_test

    return diff


### 2.2 Custom Rolling Diff function for Lasso Updating method
def rolling_diff_Lasso(d_var, ind_vars, forecast_start, wk=8, window=5, cvs=5):
    '''
    This function calculates the weekly prediction error of the Lasso Model.
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_start: current week which we base on to forecast 
        4. wk: 8 or 4 according to which vars we are interested in
        5. zero: True if 0 Specification model
        6. window: The backward looking length for coefficient updating
    Outputs:
        1. diff: difference between forecast and real observation
    '''


    ### 0. Read the dataset
    data=data_set(d_var)
    
    ### 1. Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    
    ### 2. Shift x to match the lag 
    lag_vars = ind_var_list(d_var, weeks=wk)
    # trend and WIPIyoy will not lag
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_{}wk'.format(wk))
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
    
    ### 3. Prepare LHS data for training and testing
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### 4. Prepare RHS data for training and testing   
    # 4.1 Add PCA series
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)

    # 4.2 Select the train set (first few rows) and the test set (the last row)
    data_xtrain= data_x_pca.iloc[:np.sum(date_update_range),:]
    data_xtest = data_x_pca.iloc[[-1],:]   
    
    ### 5. Lasso Update and Prediction Process
    # 5.1 First, use grid search to choose the best Penalty 
    #     Then, run lasso using that panalty regression and update, 
    #     here, sm.add_constant adds a constant to RHS
    x_train=data_xtrain.loc[:,ind_vars]
    X_train=sm.add_constant(x_train)
    ## Choose the correct LHS var according to forecasting duration
    if wk==8:
        y_train=data_ytrain[d_var+'_t8']
    else:
        y_train=data_ytrain[d_var+'_t4']
    ## Set up Lasso instance and grid search for penalty coefficient
    pre_model=Lasso()
    param_grid=[{'alpha':np.linspace(0,2,40)}]
    grid_search = GridSearchCV(pre_model, param_grid, cv=cvs, scoring='neg_mean_squared_error')
    train_xy=pd.concat([X_train,y_train],axis=1).dropna()
    y_train=train_xy.iloc[:,-1]
    X_train=train_xy.iloc[:,0:-1]
    grid_search.fit(X_train, y_train)
    best_lambda=grid_search.best_params_['alpha']
    ## Update the coefficients using the selected penalty 
    reg=Lasso(alpha=best_lambda)
    reg.fit(X_train,y_train)
    x_test=data_xtest.loc[:,ind_vars]
    X_test=sm.add_constant(x_test, has_constant='add')
    
    ## Predict and record the difference
    # Choose the correct LHS var according to forecasting duration
    if wk==8:
        y_test=data_ytest[d_var+'_t8']
    else:
        y_test=data_ytest[d_var+'_t4']
    # Prepare proper test data and predict   
    test_xy=pd.concat([X_test,y_test],axis=1).dropna()
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    
    ### 6. Return if valid o/w return np.nan
    if len(y_test)==0:
        return [np.nan]
    else:
        diff = reg.predict(X_test) - y_test

    return diff

### 2.3 Custom Rolling Diff function for Forward Model 
def rolling_diff_forward(data, d_var, ind_vars, forecast_start, wk=8, window=5, cvs=5):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_start: current week which we base on to forecast 
        4. wk: 8 or 4 according to which vars we are interested in
        5. zero: True if 0 Specification model
        6. window: The backward looking length for coefficient updating
    Outputs:
        1. diff: difference between forecast and real observation
        2. coef list: list of the coefficients for the forward selected variables
    '''
    
    ### Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    
    ### Shift x to match y and set Newey-West max lag
    lag_vars = ind_var_list(d_var, weeks=wk)
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_{}wk'.format(wk))
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
    
    ### train suffix is the update window, test suffix is the forecast window
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### Update PCA series in current lookback window   
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)

    ### Get the training data with PCA series augmented
    data_xtrain= data_x_pca.iloc[:np.sum(date_update_range),:]
    data_xtest = data_x_pca.iloc[[-1],:]   
    
    ### First, use grid search to choose the best Penalty 
    ### Then, run lasso using that panalty regression and update, here, sm.add_constant adds a constant to RHS
    x_train = data_xtrain.loc[:,ind_vars]
    X_train = sm.add_constant(x_train)
    if wk==8:
        y_train=data_ytrain[d_var+'_t8']
    else:
        y_train=data_ytrain[d_var+'_t4']
    pre_model=Lasso()
    param_grid=[{'alpha':np.linspace(0,2,40)}]
    grid_search = GridSearchCV(pre_model, param_grid, cv=cvs, scoring='neg_mean_squared_error')
    
    train_xy=pd.concat([X_train,y_train],axis=1).dropna()
    y_train=train_xy.iloc[:,-1]
    X_train=train_xy.iloc[:,0:-1]
    grid_search.fit(X_train, y_train)

    best_lambda=grid_search.best_params_['alpha']
        
    reg=Lasso(alpha=best_lambda)
    reg.fit(X_train,y_train)
    x_test=data_xtest.loc[:,ind_vars]
    X_test=sm.add_constant(x_test, has_constant='add')
    
    ### Predict and record the difference
    if wk==8:
        y_test=data_ytest[d_var+'_t8']
    else:
        y_test=data_ytest[d_var+'_t4']
        
    test_xy=pd.concat([X_test,y_test],axis=1).dropna()
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    
    ### Return the prediction error and the coefficients
    if len(y_test)==0:
        return [np.nan], reg.coef_.tolist()
    else:
        diff = reg.predict(X_test) - y_test

    return diff, reg.coef_.tolist()

### 2.4 Simple rolling diff function, save time reading dataset Return difference series between forecasted and real value
def rolling_diff(data, d_var, ind_vars, forecast_start):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_start: current week which we base on to forecast 
    Outputs:
        1. diff: difference between forecast and real observation
    '''
    
    ### Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=8, update_window=5)
    
    ### Shift x to match y and set Newey-West max lag
    lag_vars = ind_var_list(d_var,weeks=8)
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_8wk')
    data_x=data.copy()
    data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    ### train suffix is the update window, test suffix is the forecast window
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### Update PCA series in current lookback window   
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)

    data_xtrain= data_x_pca.iloc[:np.sum(date_update_range),:]
    data_xtest = data_x_pca.iloc[[-1],:]   
    
    ### First, use grid search to choos the best Penalty 
    ### Then, run lasso using that panalty regression and update, here, sm.add_constant adds a constant to RHS
    x_train=data_xtrain.loc[:,ind_vars]
    X_train=sm.add_constant(x_train)
    y_train=data_ytrain[d_var+'_t8']

    model=lm.OLS(y_train, X_train, missing='drop')
    reg=model.fit(cov_type='HAC', cov_kwds={'maxlags':5})
    x_test=data_xtest.loc[:,ind_vars]
    X_test=sm.add_constant(x_test, has_constant='add')
    
    ### Predict and record the difference
    y_test=data_ytest[d_var+'_t8']
        
    test_xy=pd.concat([X_test,y_test],axis=1).dropna()
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    
    if len(y_test)==0:
        return np.nan
    else:
        diff = y_test - reg.predict(X_test)

    return diff.values[0]

### 2.5 Make use of pre-generated coefficients for each fixced model in stable model calculation, save dozens of time | Prediction function for each week
def prediction_one_week(data, d_var, ind_vars, model_coefs, forecast_week):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_week: current week which we base on to forecast 
        4. d_var: dependent Var
        5. model_coefs: pre-generated model coefficient pd DataFrame
    Outputs:
        1. diff: difference between forecast and real observation
    '''

    ### Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_week, wk=8, update_window=5)
    
    ### Shift x to match y and set lag
    lag_vars = ind_var_list(d_var,weeks=8)
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_8wk')
    data_x=data.copy()
    data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    
    ### train suffix is the update window, test suffix is the forecast window
    data_ytest = data[date_test_range]
          
    ### Update PCA values        
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)
    
    data_xtest = data_x_pca.iloc[[-1],:]  
    
    ### Predict and record the difference
    y_test=data_ytest[d_var+'_t8']
    
    ### Get all qualified models and calculate the prediction

    cfs = model_coefs[forecast_week]
    ### Get series for each model
    x_test=data_xtest.loc[:,ind_vars.split(', ')]
    ### Add constant
    X_test=sm.add_constant(x_test, has_constant='add')      
    test_xy=pd.concat([X_test,y_test],axis=1)
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    ### Predict and save
    prediction = sum(X_test.values.squeeze()[i]*cfs[i] for i in range(3))

    ### Real Y value
    y_true = y_test.values[0]
    ### Calculate prediction error and the error of the constant model
    return [y_true-prediction], [rolling_diff(data, d_var, [], forecast_week)]

### 2.4 Custom Rolling Diff function for Stability Model (coef generating phase)
def rolling_diff_stability_coef(data, d_var, ind_vars, forecast_start, wk=8, window=5, cvs=5):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_start: current week which we base on to forecast 
        4. wk: 8 or 4 according to which vars we are interested in
        5. window: The backward looking length for coefficient updating
    Outputs:
        1. diff: difference between forecast and real observation
        2. a list with all the coefficients and R2 value of the current predicting period
    '''
    
    ### 0. Read the dataset
    data=data_set(d_var)
    
    ### 1. Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    
    ### 2. Shift x to match the lag 
    lag_vars = ind_var_list(d_var, weeks=wk)
    # trend and WIPIyoy will not lag
    lag_vars.remove('trend')
    lag_vars.remove('WIPImom_{}wk'.format(wk))
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
    
    ### 3. Prepare LHS data for training and testing
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### 4. Prepare RHS data for training and testing   
    # 4.1 Add PCA series
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca)

    # 4.2 Select the train set (first few rows) and the test set (the last row)
    data_xtrain= data_x_pca.iloc[:np.sum(date_update_range),:]
    data_xtest = data_x_pca.iloc[[-1],:]   
    
    ### 5. Lasso Update and Prediction Process
    # 5.1 First, use grid search to choose the best Penalty 
    #     Then, run lasso using that panalty regression and update, 
    #     here, sm.add_constant adds a constant to RHS
    x_train=data_xtrain.loc[:,ind_vars]
    X_train=sm.add_constant(x_train)
    ## Choose the correct LHS var according to forecasting duration
    if wk==8:
        y_train=data_ytrain[d_var+'_t8']
    else:
        y_train=data_ytrain[d_var+'_t4']
    ## Set up Lasso instance and grid search for penalty coefficient
    pre_model=Lasso()
    param_grid=[{'alpha':np.linspace(0,2,40)}]
    grid_search = GridSearchCV(pre_model, param_grid, cv=cvs, scoring='neg_mean_squared_error')
    train_xy=pd.concat([X_train,y_train],axis=1).dropna()
    y_train=train_xy.iloc[:,-1]
    X_train=train_xy.iloc[:,0:-1]
    grid_search.fit(X_train, y_train)
    best_lambda=grid_search.best_params_['alpha']
    ## Update the coefficients using the selected penalty 
    reg=Lasso(alpha=best_lambda)
    reg.fit(X_train,y_train)
    x_test=data_xtest.loc[:,ind_vars]
    X_test=sm.add_constant(x_test, has_constant='add')
    
    ## Predict and record the difference
    # Choose the correct LHS var according to forecasting duration
    if wk==8:
        y_test=data_ytest[d_var+'_t8']
    else:
        y_test=data_ytest[d_var+'_t4']
    # Prepare proper test data and predict   
    test_xy=pd.concat([X_test,y_test],axis=1).dropna()
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    
    ### 6. Return the prediction error, the coefficients and the R2
    if len(y_test)==0:
        return [np.nan], [reg.intercept_]+reg.coef_.tolist()[1:]+[reg.score(X_train,y_train)]
    else:
        diff = reg.predict(X_test) - y_test

    return diff, [reg.intercept_]+reg.coef_.tolist()[1:]+[reg.score(X_train,y_train)]