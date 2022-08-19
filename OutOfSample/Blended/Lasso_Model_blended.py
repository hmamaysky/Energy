#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""

Created on Tue Apr 28 17:36:26 2020

@author: Hongyu Wu

This file performs the OOS 1-1 and 2-2 Lasso update model.

1. The main algorithm first select a significant 1-1 or 2-2 pair (#of baseline var-#of text var)
   weekly based on a 5-yr lookback window. The siginificant variables are those with the top R-squared value
   in the univariate regression with specification : d_var = const + a*var + error.
2. After the selection phase, the coefficients of the outstanding variables will be updated
   using Lasso regression based on the same 5-yr lookback window.
3. Then forecast 4 or 8 wk ahead and record the prediction error
4. March 1 week ahead and calculate the RMSE after reaching all the available weeks

"""

# %% 0. Importing Pachages

import pandas as pd
import sys
import pickle
import concurrent.futures
from OOSfuncs import *


# %% 1. Defining Functions
### See OOSfuncs.py for all the functions ###

# %% 2. Main Function
def main(d_var):
    """
    This main function carry out all the steps of the Lasso Model algorithm:
        1. At each window, look backward for 5 years and select RHS variables.
        2. Update coefficients of the selected ones using that 5 year window again with Lasso regression.
        3. Forecast next month or two month (this can be controlled with the sys.arg[1] or the var 'weeks')
        4. Calculate discrepancy between prediction and observation
        5. Calculate and save RMSE for each model and each LHS variables
    Note:
        1. A DataFrame that is 4*8 (4 models, 8 LHS vars)
        2. 4 models : constant, 1-1, 1-0, and 0-1
    """
    ### 0. Set Parameters
    # Read Sys Args
    weeks=int(sys.argv[1])
    frequency=int(sys.argv[2])
    no_varibles=int(sys.argv[3])
    cv=int(sys.argv[4])
    # weeks=8
    # frequency=1
    # no_varibles=2
    # cv=5    
    # d_var list
    d_vars = [d_var]

    # Lookback window (/yr)
    updating_windows = [5]

    ### 1. Main Algorithm
    for updating_window in updating_windows:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print("Update Window: ", updating_window)
        
        # 1.1 Result holder: 
        # a dictionary with each dependent variable as key,
        # and a list of rmses for each model as value
        rmse_result={}
        rmse_series={}
        prediction={}
        for d_var in d_vars:
            print('-------------------')
            print(d_var)
            
            # 1.2 More result holder for the prediction differences
            #     Model names are suggested by the var name
            constant_diff_list =[]
            base_diff_list =[]
            text_diff_list =[]
            full_diff_list =[]
            blended_diff_list =[]

            constant_pred_list =[]
            base_pred_list =[]
            text_pred_list =[]
            full_pred_list =[]
            blended_pred_list =[]
            
            #     List for cumulative rmse
            constant_rmse_list =[]
            base_rmse_list =[]
            text_rmse_list =[]
            full_rmse_list =[]
            blended_rmse_list =[]
                        
            # 1.3 Get the exact date of each window (monthly)
            #     Here, we ensure everytime we train, there are enough observations for backward looking
            time_col = data_set(d_var)['date'] 
            time_lower=time_col[0]
            # time_lower= pd.Timestamp('2014-05-03')
            test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
            
            # 1.4 Main algorithm
            # First, at each window, train and select variables
            # Second, update coefficients on those variables
            #   (Note that above two steps are done using 5-year backward looking window)
            # Then, forecast and get its difference from real observations
            # Finally, record them in a list
            for week in test_week_list:
                print(week)
                
                # Train and select RHS variables
                full_vars, base_vars, text_vars = select_significant(d_var, week, wk=weeks, ns=no_varibles)
                
                # Calculate MSE and RMSE
                constant_diff, constant_pred = rolling_diff_Lasso(d_var, [], week, wk=weeks, window=updating_window, cvs=cv)
                base_diff, base_pred = rolling_diff_Lasso(d_var, base_vars, week, wk=weeks, window=updating_window, cvs=cv)
                text_diff, text_pred = rolling_diff_Lasso(d_var, text_vars, week, wk=weeks, window=updating_window, cvs=cv)
                full_diff, full_pred = rolling_diff_Lasso(d_var, full_vars, week, wk=weeks, window=updating_window, cvs=cv)
                w=0 # set weights
                blended_diff = w*constant_diff + (1-w)*full_diff
                blended_pred = w*constant_pred + (1-w)*full_pred
               
                # Append results to lists              
                constant_diff_list.extend(constant_diff)                
                base_diff_list.extend(base_diff)    
                text_diff_list.extend(text_diff) 
                full_diff_list.extend(full_diff)
                blended_diff_list.extend(blended_diff)
                
                constant_pred_list.extend(constant_pred)                
                base_pred_list.extend(base_pred)    
                text_pred_list.extend(text_pred) 
                full_pred_list.extend(full_pred)
                blended_pred_list.extend(blended_pred)

                constant_rmse_list.append(RMSE(constant_diff_list))
                base_rmse_list.append(RMSE(base_diff_list))
                text_rmse_list.append(RMSE(text_diff_list))
                full_rmse_list.append(RMSE(full_diff_list))           
                blended_rmse_list.append(RMSE(blended_diff_list))
                
            # 1.5 Store RMSE ratios of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(base_diff_list)/RMSE(constant_diff_list), RMSE(text_diff_list)/RMSE(constant_diff_list),
                          RMSE(full_diff_list)/RMSE(constant_diff_list), RMSE(text_diff_list)/RMSE(base_diff_list),
                          RMSE(full_diff_list)/RMSE(base_diff_list)]
            rmse_series[d_var] = {'const_rmse':constant_rmse_list,
                                  'base_rmse':base_rmse_list,
                                  'text_rmse':text_rmse_list,
                                  'full_rmse':full_rmse_list,
                                  'blended_rmse':blended_rmse_list}
            
            prediction[d_var] = {'const_pred':constant_pred_list,
                                  'base_pred':base_pred_list,
                                  'text_pred':text_pred_list,
                                  'full_pred':full_pred_list,
                                  'blended_pred':blended_pred_list}
        # 2. Save and Return 
        # Create a DataFrame to save results, more convenient for further interpretation
        rmse_df = pd.DataFrame(rmse_result)
        rmse_df.index=['Base Model', 'Text Model', 'Full Model', 'Text Model (on Base)', 'Full Model (on Base)']

    return rmse_df, rmse_series, prediction

# %% 3. Main Process
if __name__ == '__main__':
    
    ### 3.1 System Args
    forecasting_week=int(sys.argv[1])
    update_frequency=int(sys.argv[2])
    no_variables=int(sys.argv[3])
    cv_fold=int(sys.argv[4])
    # forecasting_week=8
    # update_frequency=1
    # no_variables=1
    # cv_fold=10
    ### 3.2 File Naming Strings
    if    forecasting_week == 4: file_suffix = '4wk'
    else: file_suffix = '8wk'
        
    if    update_frequency == 1: file_start = 'weekly'
    else: file_start = 'monthly'
        
    if    no_variables  == 1: file = 'oneandone'
    elif  no_variables == 2: file = 'twoandtwo'
    else: file = 'threeandthree'

    ### 3.3 Use concurrent calculation for the main algorithm
    # result holders for each process
    rmse=pd.DataFrame()
    rmse_seriess=dict()
    prediction_series=dict()
    # 8 dependent variables
    d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    # Main algorithm and save the results
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(main,d_vars)
        for result in results:
            rmse=pd.concat([rmse,result[0]], axis=1)
            rmse_seriess.update(result[1])
            prediction_series.update(result[2])
            
    ### 3.4 Save the results to proper directory
    rmse.to_excel('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                    +file+'/'+file_start+'/Lasso_'+str(cv_fold)+'fold_'+file_suffix+'_t.xlsx')
    pickle.dump(rmse_seriess,open('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                    +file+'/'+file_start+'/Lasso_'+str(cv_fold)+'fold_'+file_suffix+'_t_rmse.p','wb'))
    pickle.dump(prediction_series,open('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                    +file+'/'+file_start+'/Lasso_prediction_'+str(cv_fold)+'fold_'+file_suffix+'_t.p','wb'))
