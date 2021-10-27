#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 21:54:04 2020

@author: Hongyu Wu

This file performs the OOS 1-1 and 2-2 OLS update model.

1. The main algorithm first select a significant 1-1 or 2-2 pair (#of baseline var-#of text var)
   weekly based on a 5-yr lookback window. The siginificant variables are those with the top R-squared value
   in the univariate regression with specification : d_var = const + a*var + error.
2. After the selection phase, the coefficients of the outstanding variables will be updated
   using OLS regression based on the same 5-yr lookback window.
3. Then forecast 4 or 8 wk ahead and record the prediction error
4. March 1 week ahead and calculate the RMSE after reaching all the available weeks
5. The outputs include the RMSEs of each model and the time series of the selected significant variables
"""

# %% 0. Importing Pachages

import pandas as pd
import sys
from OOSfuncs import *

# %% 1. Defining Functions
### See OOSfuncs.py for all the functions ###

# %% 2. Main Process
def main():
    """
    This main function carry out all the steps of the 1-and-1 algorithm:
        1. At each window, look backward for 5 years and select RHS variables.
        2. Update coefficients of the selected ones using that 5 year window again.
        3. Forecast next month or two month (this can be controlled in the "get_test_row" function)
        4. Calculate discrepancy between prediction and observation
        5. Calculate and save RMSE for each model and each LHS variables
    Note:
        1. The output is a list of dataFrames, each dataFrame are 4*8 (4 models, 8 LHS vars)
        2. The length of output list is in accord with the list "updating_windows", which can be adjusted to test the impact of different look back periods
    """
    ### 0. Set Parameters
    # Read Sys Args
    weeks=int(sys.argv[1])
    frequency=int(sys.argv[2])
    no_varibles=int(sys.argv[3])
    # weeks=8
    # frequency=1
    # no_varibles=1  
    # d_var list
    d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    # d_vars = ['DInv', 'DProd']
    # Result holder and the lookback window (/yr)
    final_rmse=[]
    updating_windows = [5]
    
    ### 1. Main Algorithm
    for updating_window in updating_windows:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print("Update Window: ", updating_window)
         
        # 1.1 Result holder: 
        # a dictionary with each dependent variable as key,
        # and a list of rmses for each model as value.
        # other two pd.DF are for the time series of the selected variables
        rmse_result={}
        text_var_df=pd.DataFrame()
        base_var_df=pd.DataFrame()
        for d_var in d_vars:
            print('-------------------')
            print(d_var)
            
            # 1.2 More result holder for the prediction differences
            #     Model names are suggested by the var name
            constant_diff_list =[]
            base_diff_list =[]
            text_diff_list =[]
            full_diff_list =[]
            
            #    var list by week, to be saved in the pd.DF created before
            base_var_list={}
            text_var_list={}
            
            # 1.3 Get the exact date of each window (monthly)
            #     Here, we ensure everytime we train, there are enough observations for backward looking
            time_col = data_set(d_var)['date'] 
            test_week_list = [time for time in time_col if time>=time_col[0]+pd.Timedelta(str(7*updating_window*52)+'days')][::frequency]
            
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
                constant_diff = rolling_diff_OLS(d_var, [], week, wk=weeks, window=updating_window)
                base_diff = rolling_diff_OLS(d_var, base_vars, week, wk=weeks, window=updating_window)
                text_diff = rolling_diff_OLS(d_var, text_vars, week, wk=weeks, window=updating_window)
                full_diff = rolling_diff_OLS(d_var, full_vars, week, wk=weeks, window=updating_window)
                
                
                # Append results to lists              
                constant_diff_list.extend(constant_diff)                
                base_diff_list.extend(base_diff) 
                text_diff_list.extend(text_diff) 
                full_diff_list.extend(full_diff)
                
                # Save var selected according to the # of selected vars
                if no_varibles ==1:
                    base_var_list[str(week.to_pydatetime())[0:10]]=base_vars[0]
                    text_var_list[str(week.to_pydatetime())[0:10]]=text_vars[0]
                elif no_varibles ==2:
                    base_var_list[str(week.to_pydatetime())[0:10]]='{},{}'.format(base_vars[0], base_vars[1])
                    text_var_list[str(week.to_pydatetime())[0:10]]='{},{}'.format(text_vars[0],text_vars[1])
                else:
                    base_var_list[str(week.to_pydatetime())[0:10]]='{},{},{}'.format(base_vars[0], base_vars[1], base_vars[2])
                    text_var_list[str(week.to_pydatetime())[0:10]]='{},{},{}'.format(text_vars[0],text_vars[1], text_vars[2])
                    
            # 1.5 Store RMSEs and the selected vars of different models for each LHS variable  
            rmse_result[d_var]=[RMSE(constant_diff_list), RMSE(base_diff_list),
                         RMSE(text_diff_list), RMSE(full_diff_list)]
            text_var_df=pd.concat([text_var_df,pd.DataFrame(text_var_list, index=[d_var])])
            base_var_df=pd.concat([base_var_df,pd.DataFrame(base_var_list, index=[d_var])])
 
        # 2. Save and Return 
        # Create a DataFrame to save results, more convenient for further interpretation
        rmse_df = pd.DataFrame(rmse_result)
        rmse_df.index=['Constant', 'Base Model', 'Text Model', 'Full Model']
        final_rmse.append(rmse_df)

    return final_rmse, text_var_df, base_var_df

# %% 3. Main Process
if __name__ == '__main__':
    ### 3.1 System Args
    forecasting_week=int(sys.argv[1])
    update_frequency=int(sys.argv[2])
    no_variables=int(sys.argv[3])
    # forecasting_week=8
    # update_frequency=1
    # no_variables=1    
    ### 3.2 File Naming Strings
    if    forecasting_week == 4: file_suffix = '4wk'
    else: file_suffix = '8wk'
        
    if    update_frequency == 1: file_start = 'weekly'
    else: file_start = 'monthly'
        
    if    no_variables  == 1: file = 'oneandone'
    elif  no_variables == 2:file = 'twoandtwo'
    else: file = 'threeandthree'
    ### 3.3 Perform Main algorithm
    rmse, text_var_df, base_var_df=main()
    
    ### 3.4 Save the results to proper directory
    rmse[0].to_excel('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                +file+'/'+file_start+'/'+file_start+'_'+file_suffix+'.xlsx')
    text_var_df.to_csv('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                       +file+'/'+file_start+'/text_vars_'+file_start+'_'+file_suffix+'.csv')
    base_var_df.to_csv('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/parsimonious/'
                       +file+'/'+file_start+'/base_vars_'+file_start+'_'+file_suffix+'.csv')