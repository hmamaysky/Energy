#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 14:57:25 2020

@author: Hongyu Wu

Description:
------------
    This new model takes average predictions of the selected models as an estimate of future return.
    The selected models are those with stable coefficients within past few (parameter) years.
    Stable means that none of the coefficients hit zero or switched signs during the inspecting window.
    If there is not such a model, skip the time stamp and march to the next week.
    Note that the coefficients are calculated beforehead to save the running time of this process.

Output: RMSE ratios of the time varying model on the corresponding constant models.
        Proportion of weeks with stable candidates.
        Average proportion of stable candidates weekly.
        
Note:
-----
    Please change the directories in the main function and the main process to proper ones before running the program.

"""

# %% 0. Importing Packages
import pandas as pd
import pickle
import sys
import os
import heapq
import functools
import concurrent.futures
import numpy as np
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
from OOSfuncs import *
# %% 1. Defining Functions
### See OOSfuncs.py for all the other functions not defined here ###

### Average function with memory
def average():
    """return the average of all the numbers entered up to current call"""
    all_elements=[]
    def take_ave(num):
        all_elements.append(num)
        summing = sum(all_elements)
        number = len(all_elements)
        return summing/number
    return take_ave

### Return difference series between forecasted and real value
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

### Prediction function for each week
def prediction_one_week(data, d_var, forecast_week, selected_models, constant_fillin=False):
    '''
    Inputs:
        1. d_var: dependent variable
        2. ind_vars: prescribed model
        3. forecaset_week: current week which we base on to forecast 
        4. selected_models:
        5. constant_fill: 
    Outputs:
        1. diff: difference between forecast and real observation
    '''
    ### 0. Check if there is any candidate model for the current forecasting week
    ###    if not, check if constant_fillin is True, if still not, return nan and nan
    if any(selected_models):
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
        predictions = []
        for model, coefs in selected_models.items():
            ### Get series for each model
            x_test=data_xtest.loc[:,model.split(', ')]
            ### Add constant
            X_test=sm.add_constant(x_test, has_constant='add')      
            test_xy=pd.concat([X_test,y_test],axis=1)
            y_test=test_xy.iloc[:,-1]
            X_test=test_xy.iloc[:,0:-1]
            ### Predict and save
            prediction = sum(X_test.values.squeeze()[i]*coefs[i] for i in range(3))
            predictions.append(prediction)
        
        ### Get rid of the na predictions (results from missing value in either x or y)
        predictions = pd.Series(predictions).dropna()
        
        ### If no prediction, and constant_fillin is True, return constant model error
        ### else, return nan, nan
        if len(predictions)==0:
            if constant_fillin:
                const_diff = rolling_diff(data, d_var, [], forecast_week)
                return const_diff, const_diff
            else:
                return np.nan, np.nan
        
        ### Similar check as above, but one predicted variable
        if len(y_test)==0:
            if constant_fillin:
                const_diff = rolling_diff(data, d_var, [], forecast_week)
                return const_diff, const_diff
            else:
                return np.nan, np.nan
        ### Real Y value
        y_true = y_test.values[0]
        ### Calculate prediction error and the error of the constant model
        return y_true-np.mean(predictions), rolling_diff(data, d_var, [], forecast_week)
    
    elif constant_fillin:
        ### return const model diff if no valid stable candidates and if constant_fillin is set to be True
        const_diff = rolling_diff(data, d_var, [], forecast_week)
        return const_diff, const_diff
    
    else:
        ### return nan, nan if no stable candidates and constant_fillin is not required
        return np.nan, np.nan
        
# %% 2. Main function
def main(d_var, lookback=3):
    # Read the whole dataset
    data = data_set(d_var)
    
    # Change dir to read in selected models and coefs
    wkdir = '/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/processed'
    os.chdir(wkdir)
    
    # Read in all the selected models and coefs
    selected_models=pickle.load(open(str(lookback)+'yrback_selected_models_with_coefs.p','rb'))
    selected_models = selected_models[d_var]
    total_forecast=0; actual_forecast=0; average_model=0
    ave=average()
    for _,j in selected_models.items():
        if any(j):average_model=ave(len(j))    
    
    # Get all the forecast weeks
    time_col = data_set(d_var)['date']
    time_lower = time_col[0]
    week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*(5+lookback)*52)+'days')][::1]

    
    # Create list to hold raw results
    timevarying_results=[]
    const_sub_results=[]
    constfillin_results=[]
    constant_results=[]
    
    # Predict and calculate the error for the varying model and the constant model
    for forecast_week in week_list:
        timevarying_error, constant_error_sub = prediction_one_week(data, d_var, forecast_week, selected_models[forecast_week])
        constfillin_error, _ = prediction_one_week(data, d_var, forecast_week, selected_models[forecast_week], constant_fillin=True)
        constant_error = rolling_diff(data, d_var, [], forecast_week)
        # save the prediction error of the const_fill model, plain model and the constant model of each respectively
        # some mismatch on the prediction horizon might occur in the constant models of the two specifications
        # please refer to the .doc file in the same directory for details
        const_sub_results.append(constant_error_sub)
        timevarying_results.append(timevarying_error)
        constfillin_results.append(constfillin_error)
        constant_results.append(constant_error)
    # Calculate related stats
    actual_forecast = len(pd.Series(timevarying_results).dropna())
    total_forecast = len(pd.Series(constant_results).dropna())
    # Save the results
    results_df = pd.DataFrame({d_var:[RMSE(timevarying_results)/RMSE(const_sub_results), RMSE(timevarying_results)/RMSE(constant_results), RMSE(constfillin_results)/RMSE(constant_results)
                                , '{}/{}'.format(actual_forecast, total_forecast), average_model/741]},
                               index=['{:01d}yr Stable Model sub'.format(lookback), '{:01d}yr Stable Model all'.format(lookback), 'Const Fillin Model', 
                                      'Stable/Constant Pred Weeks', 'Model Proportion'])
    return results_df

# %% 3. Main process
if __name__=='__main__':
    # Get lookback window from sys arg
    lookback = int(sys.argv[1])
#    lookback=3
    # Kick off the calculation by parallelization
    dvar_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    rmse_results = pd.DataFrame()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(functools.partial(main, lookback=lookback), dvar_list)
        for result in results:
            rmse_results = pd.concat([rmse_results, result], axis=1)
    # Save the results
    rmse_results.to_excel('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/outputs/'+
                          str(lookback)+'yrstable_model.xlsx')
    
    