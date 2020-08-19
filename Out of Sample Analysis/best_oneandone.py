#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:03:22 2020

This code pick the best one and one model for each dependent variable based on RMSE.

The RMSE calculation is conducted using the other code, this is just the summary part.

This model fixed the 1-1 pair between text and baseline variables at first,
 and update the coefficients using 5-yr backward looking window as well as making prediction weekly.
 The results are rather promising because we are optimizing the RMSE using all the data.

@author: Hongyu Wu
"""
# %% 0. Importing Packages
import pandas as pd
import numpy as np
import os

# %% 1. Defining Functions

## 1.1 Concatenating all the results in the fixed model implementation.
def data_concat():
    ''' This function concat all the results. '''
    ## Take the last file as the starting rows
    ## The Constant row is not needed for the other ones because they are identical
    ## Special cases involving sdf and ovx vars wil be handled later
    data = pd.read_excel(files[-1])
    for i in range(len(files)-1):
        temp = pd.read_excel(files[i])
        data = pd.concat([data, temp.iloc[1:,:]])
    data.rename(columns={'Unnamed: 0': 'model'}, inplace = True)
    data=data[data['model']!='Constant']
    return data

## 1.2 Extract all the models that beats their corresponding constant model.
def less_than_const():
    '''
    This function output the models that beats the constant model and their MSE ratios against the constant.
    There are 8 dependent vars, so 8 .csv files are produced.
    For ovx, and sdf series, the time horizon is different from the regular variables, so the calculation,
    and the comparison should be separately treated.
    '''
    # Save the constant model rmse into a dictionary for different cases
    constmodel_RMSEs=dict()
    constmodel_RMSEs['ovx'] = pd.read_excel('ovx_cl1_Lasso_10fold_8wk_M.xlsx').loc[0,:]
    constmodel_RMSEs['sdf'] = pd.read_excel('RPsdf_growing_Lasso_10fold_8wk_M.xlsx').loc[0,:]
    constmodel_RMSEs['regular'] = pd.read_excel('trend_Lasso_10fold_8wk_M.xlsx').loc[0,:]
    # Read all the data
    data = data_concat()
    # Separate data into ovx related, sdf related and regular models
    data_ovx = data[data['model'].str.contains('ovx')]
    data_sdf = data[data['model'].str.contains('sdf')& ~data['model'].isin(data_ovx['model'].tolist())]
    data_reg = data[~(data['model']).isin(data_ovx['model'].tolist()+data_sdf['model'].tolist())]
    # Start the main selection and calculation
    var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot',
       'DOilVol', 'DInv', 'DProd']
    for var in var_list:
        # For each cases (ovx, sdf, regular), select the model that beats the constant model
        # and calculate the MSE ratio against the constant model
        temp_ovx = data_ovx[data_ovx[var]<constmodel_RMSEs['ovx'][var]].loc[:,['model',var]]
        temp_ovx['ratio'] = (temp_ovx.loc[:,var]/constmodel_RMSEs['ovx'][var])**2
        temp_sdf = data_sdf[data_sdf[var]<constmodel_RMSEs['sdf'][var]].loc[:,['model',var]]
        temp_sdf['ratio'] = (temp_sdf.loc[:,var]/constmodel_RMSEs['sdf'][var])**2
        temp_reg = data_reg[data_reg[var]<constmodel_RMSEs['regular'][var]].loc[:,['model',var]]
        temp_reg['ratio'] = (temp_reg.loc[:,var]/constmodel_RMSEs['regular'][var])**2
        # integrate all the results
        temp = pd.concat([temp_ovx, temp_sdf, temp_reg], axis=0)
        # drop duplicates by utilizing frozenset
        temp['set'] = temp.loc[:,'model'].apply(lambda x: frozenset(x))
        temp = temp.drop_duplicates(subset='set')
        temp = temp.drop('set', axis=1)
        # Sort the df by MSE ratio
        temp.sort_values(by='ratio', ascending=True, inplace=True)
        # output the results
        temp.to_csv('../results/'+var+'.csv',index=False)
        # output the results for Yearly Diff calculation
        temp.loc[:,['model','ratio']].to_csv('../summary/'+var+'.csv',index=False)

# %% 2. Executing the whole process
if __name__ == "__main__":
    ## Set the proper directory (the one containing all the raw results of the Part 1-3 of the fixed model implementation)
    wkdir = 'Please set your results'
    os.chdir(wkdir)
    ## Get all the files
    files = os.listdir()
#    ## remove irrelevant files if necessary
#    remove_list = ['results', 'code', 'prediction_real',
#                   'pred_and_real.p', 'pred_and_real_new.p',
#                   'summary', '.DS_Store']
#    for i in remove_list:
#        try:     files.remove(i)
#        except:  pass
    ## The "main" process
    less_than_const()
    