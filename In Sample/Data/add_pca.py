#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:29:14 2020

@author: billwu
"""

import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA


def add_PCA(df, sent_vars, freq_vars, all_vars):
    data_temp = df.copy()
    pca = PCA(n_components=1)
    
    ## PCA columns
    data_temp['PCAsent'] = np.nan
    data_temp['PCAfreq'] = np.nan
    data_temp['PCAall'] = np.nan
    
    ## The data series starts from the 5th row
    data_temp.loc[4:, 'PCAsent'] = pca.fit_transform(data_temp[sent_vars][4:])
    data_temp.loc[4:, 'PCAfreq'] = pca.fit_transform(data_temp[freq_vars][4:])
    data_temp.loc[4:, 'PCAall'] = pca.fit_transform(data_temp[all_vars][4:])

    return data_temp

def main():
    
    ## Set the wkdir
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/'
    wkdir += 'Energy Project/Analysis/data/transformed_data'
    os.chdir(wkdir)
    
    ## Read the data
    data_price = pd.read_stata('transformed_data_prices_v14.dta')
    data_physical = pd.read_stata('transformed_data_physical_v14.dta')
    
    ## Specify text, sent and freq variables
    all_vars = ['artcount', 'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp']  
    all_vars_drop = all_vars.copy()
    all_vars_drop.remove('fCo')
    freq_vars = ['fCo', 'fGom', 'fEnv', 'fEpg', 'fBbl', 'fRpc', 'fEp']
    freq_vars_drop = freq_vars.copy()
    freq_vars_drop.remove('fCo')
    sent_vars = ['sCo', 'sGom', 'sEnv', 'sEpg', 'sBbl', 'sRpc', 'sEp']
    
    ## Add the first principal components into the data set
    data_price_pca = add_PCA(data_price, sent_vars, freq_vars_drop, all_vars_drop)
    data_physical_pca = add_PCA(data_physical, sent_vars, freq_vars_drop, all_vars_drop)
    
    ## Save the data set
    data_price_pca.to_stata('transformed_data_prices_v14.dta')
    data_physical_pca.to_stata('transformed_data_physical_v14.dta')
    
if __name__ == '__main__':
    main()
