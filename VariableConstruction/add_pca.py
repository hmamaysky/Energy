"""
This script creates the following variables: 

1. PCAsent
2. PCAfreq
3. PCAall

Outputs: 
1. transformed_data_physical_pca_v19.dta
2. transformed_data_prices_pca_v19.dta
"""

import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA


def add_PCA(df, df_temp, sent_vars, freq_vars, all_vars):
    data_temp = df.copy()
    pca = PCA(n_components=1)
    
    ## PCA columns
    data_temp['PCAsent'] = np.nan
    data_temp['PCAfreq'] = np.nan
    data_temp['PCAall'] = np.nan
    
    ## The data series starts from the 5th row
    data_temp.loc[3:, 'PCAsent'] = pca.fit_transform(df_temp[sent_vars][3:])
    data_temp.loc[3:, 'PCAfreq'] = pca.fit_transform(df_temp[freq_vars][3:])
    data_temp.loc[3:, 'PCAall'] = pca.fit_transform(df_temp[all_vars][3:])

    return data_temp

def main():
    
    ## Set the wkdir
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/'
    wkdir += 'Energy Project/Analysis/data/transformed_data'
    os.chdir(wkdir)
    
    ## Read the data
    data_price = pd.read_stata('transformed_data_prices_v19.dta')
    date_cols_price = [x for x in list(data_price.columns.values) if 'date' in x]
    data_price_temp = data_price.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(data_price.columns.values) if x not in date_cols_price})
    data_physical = pd.read_stata('transformed_data_physical_v19.dta')
    date_cols_physical = [x for x in list(data_physical.columns.values) if 'date' in x]
    data_physical_temp = data_physical.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(data_physical.columns.values) if x not in date_cols_physical})
    
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
    data_price_pca = add_PCA(data_price, data_price_temp, sent_vars, freq_vars_drop, all_vars_drop).rename(columns={'PCAsent':'PCAsent_Fri',
                                                                                                            'PCAfreq':'PCAfreq_Fri',
                                                                                                            'PCAall':'PCAall_Fri'})
    data_physical_pca = add_PCA(data_physical, data_physical_temp, sent_vars, freq_vars_drop, all_vars_drop).rename(columns={'PCAsent':'PCAsent_Tue',
                                                                                                            'PCAfreq':'PCAfreq_Tue',
                                                                                                            'PCAall':'PCAall_Tue'})
    
    ## Save the data set
    data_price_pca = data_price_pca.loc[:,date_cols_price+['PCAsent_Fri','PCAfreq_Fri','PCAall_Fri']]
    for dates in date_cols_price:
        data_price_pca[dates] = data_price_pca[dates].apply(lambda x:str(x)[:10])
    data_physical_pca = data_physical_pca.loc[:,date_cols_physical+['PCAsent_Tue','PCAfreq_Tue','PCAall_Tue']]
    for dates in date_cols_physical:
        data_physical_pca[dates] = data_physical_pca[dates].apply(lambda x:str(x)[:10])
#    return data_price_pca
    data_price_pca.to_stata('transformed_data_prices_pca_v19.dta',write_index=False)
    data_physical_pca.to_stata('transformed_data_physical_pca_v19.dta',write_index=False)
    
        
if __name__ == '__main__':
    main()
