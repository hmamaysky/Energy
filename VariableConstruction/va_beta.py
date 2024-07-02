"""
This script creates the following variables: 

1. cpiyr_betas - renamed as InflaBeta in the final dataset
2. dxy_betas - renamed as DolBeta in the final dataset

"""

# %%
import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
__data_dir__={'raw':'/Users/Economist/Documents/Energy Project/data/raw/',
              'result':'/Users/Economist/Documents/Energy Project/data/results/'}

# %%
def calc_beta(df_, y_var, x_var, window=60):
    ### Get a deep copy of original data set
    df = df_.copy()
    ### Get the rolling diff for x_var (CPIyoy); For DXY, we already have log-diff variable.
    if x_var == 'cpiyr': 
        df.loc[:,[x_var]] = df.loc[:,[x_var]].diff()
    
    df = df.dropna(subset=[x_var])
    ### Rolling OLS Regressions
    betas = []
    for i in range(df.shape[0]-window+1):
        temp = df.loc[i:i+60, [x_var,y_var]]
        X = sm.add_constant(temp.loc[:,[x_var]])
        Y = temp.loc[:,[y_var]]
        ols = sm.OLS(Y,X)
        ols_res = ols.fit()
        beta = ols_res.params[x_var]
        betas.append(beta)
    ### Get Months
    months = list(df['monthly'])[window-1:]  
    ### Return
    return pd.DataFrame({'monthly':months,x_var+'_betas':betas})

# %%
if __name__ == '__main__':
    ### Read the data for futures return, CPIyoy and DXY
    return_data = pd.read_stata(__data_dir__['raw']+'cpidxyfutret_v2.dta')
    ### Run the rolling regressions for Betas
    InflaBeta = calc_beta(return_data, 'FutRet1', 'cpiyr', window=60)
    DolBeta = calc_beta(return_data, 'FutRet1', 'dxy', window=60)
    ### Save it to STATA
    InflaBeta.to_stata(__data_dir__['result']+'cpi_beta_v2.dta')
    DolBeta.to_stata(__data_dir__['result']+'dxy_beta_v2.dta')