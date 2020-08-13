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
import os
import sys
import pickle
import warnings
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA
warnings.filterwarnings('ignore')

# %% 1. Defining Functions

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
def PCA_augment(data, textual_vars_drop, freq_vars_drop, sent_vars):
    '''
    This function augments the PCA series for a given lookback window
    '''
    pca = PCA(n_components=1)
    data_temp = data.copy()
    data_temp['PCAsent'] = pca.fit_transform(data_temp[sent_vars])
    data_temp['PCAall'] = pca.fit_transform(data_temp[textual_vars_drop])
    data_temp['PCAfreq'] = pca.fit_transform(data_temp[freq_vars_drop])
    return data_temp

### 1.6 Return difference series between forecasted and real value ###
def rolling_diff(data, d_var, ind_vars, forecast_start, wk=8, window=5, cvs=5):
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
    
    ### Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)
    
    ### Shift x to match y and set Newey-West max lag
    lag_vars = ind_var_list(d_var, weeks=wk)
    lag_vars.remove('trend')
    lag_vars.remove('WIPIyoy')
    data_x=data.copy()
    if wk == 8:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(8)

    else:
        data_x.loc[:,lag_vars]=data_x.loc[:,lag_vars].shift(4)
    
    ### train suffix is the update window, test suffix is the forecast window
    data_ytrain= data[date_update_range]
    data_ytest = data[date_test_range]
    
    ### Update PCA series in current lookback window   
    ## Drop fCo in textual var list and freq var list to avoid collinearity
    textual_vars_drop = ['artcount', 'entropy', 'sCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']  
    freq_vars_drop = ['fGom', 'fEnv', 'fEpg', 'fBbl', 'fRpc', 'fEp']
    sent_vars = ['sCo', 'sGom', 'sEnv', 'sEpg', 'sBbl', 'sRpc', 'sEp']
    data_x_pca = data_x[date_pca_range]
    data_x_pca = PCA_augment(data_x_pca, textual_vars_drop, freq_vars_drop, sent_vars)

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
        y_train=data_ytrain[d_var]
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
        y_test=data_ytest[d_var]
        
    test_xy=pd.concat([X_test,y_test],axis=1).dropna()
    y_test=test_xy.iloc[:,-1]
    X_test=test_xy.iloc[:,0:-1]
    
    if len(y_test)==0:
        return [np.nan], reg.coef_.tolist()
    else:
        diff = reg.predict(X_test) - y_test

    return diff, reg.coef_.tolist()

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
            constant_diff,_ = rolling_diff(data, d_var, [], week, wk=weeks, window=updating_window, cvs=cv)
            forward_diff, var_coeff = rolling_diff(data, d_var, forward_vars, week, wk=weeks, window=updating_window, cvs=cv)
            
            
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
        # Main algorith and save the results
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