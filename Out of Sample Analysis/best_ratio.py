#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 17:58:48 2020

Calculate the MSE ratio between the winning models and the constant model.
Save the results in the /summary folder, will be further used in the matrix plot process.

@author: Hongyu Wu
"""
# %% 0. Importing Packages
import pandas as pd
import numpy as np
import os

# %% 1. Defining Functions

## 1.1 Function for the mse ratio of the winning models against the constant model
def mse_ratio(cols, const_val, d_var):
    '''
    Return the mse ratio of the model against the constant model.
    Data series for Baseline vars in the special list starts later than other series, 
    so they are separately treated.
    '''
    special_list = ['ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
    base_var = cols['model'].split(',')[0]
    if base_var in special_list:
        temp = pd.read_excel('../'+base_var+'_Lasso_10fold_8wk_M.xlsx')
        return np.power(cols['rmse']/temp.loc[0, d_var],2)
    else:
        temp = pd.read_excel('../FutRet_Lasso_10fold_8wk_M.xlsx')
        return np.power(cols['rmse']/temp.loc[0, d_var],2)
    
# %% 2. Main function
def main():
    ## LHS var list
    d_vars = ['bpRet', 'DInv', 'DOilVol', 'DProd', 'DSpot', 'FutRet', 'rdsaRet', 'xomRet']
    for d_var in d_vars:
        ## Get all the winning models
        df = pd.read_csv(d_var+'.csv', header=None).rename(columns={0:'model',1:'rmse'})
        ## Get the constant model RMSE
        const_val = df.iloc[0,1]
        ## Calculate the MSE ratios
        df['ratio']=df.apply(mse_ratio, axis=1, args=(const_val, d_var))
        df = df.drop('rmse', axis=1)
        ## Sort by ratio ascendingly
        df.sort_values(by=['ratio'], ascending=True, inplace=True)
        ## Save the result
        df.to_csv('../summary/'+d_var+'.csv', header=False, index=False)

# %% 3. Executing the main process      
if __name__ == "__main__":
    ## Set the Proper Directory (where all the results from the 'best_oneandone.py' is saved)
    wkdir = 'Please Set the Proper directory'
    os.chdir(wkdir)
    ## Main Process
    main()