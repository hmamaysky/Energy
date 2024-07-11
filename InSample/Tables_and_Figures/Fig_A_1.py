#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 12:04:42 2020

@author: billwu
"""

import numpy as np
import pandas as pd
import os
import re
import warnings
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy import signal
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
warnings.filterwarnings('ignore')
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['axes.labelweight'] = 'bold'


def plot_corr(df,size=10):
    '''Plot a graphical correlation matrix for a dataframe.

    Input:
        df: pandas DataFrame
        size: vertical and horizontal size of the plot'''
    # Compute the correlation matrix for the received dataframe
    df = df.rename(columns={'BEME':'BE/ME','basismom':'BasMom','mom_fut1':'Mom','dxy_betas':'DolBeta',
                            'cpiyr_betas':'InflaBeta','hp':'HedgPres','liquidity':'liquidity',
                            'oi':'OpenInt','WIPImom_8wk':'WIPI','ovx_cl1':'ovx_diff','vix_spx':'vix_diff',
                            'StikIdx':'StkIdx'})
    corr = df.corr()
    
    # Plot the correlation matrix
    fig, ax = plt.subplots(figsize=(size, size), dpi=200)
    cax = ax.matshow(corr, cmap='RdYlGn')

    plt.xticks(range(len(corr.columns)), [re.sub('_8wk','',x) for x in corr.columns.values.tolist()], rotation=90); # get rid of the subscripts in var names
    plt.yticks(range(len(corr.columns)), [re.sub('_8wk','',x) for x in corr.columns.values.tolist()]); # get rid of the subscripts in var names
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom, top)
    
    # Add the colorbar legend
    cbar = fig.colorbar(cax, ticks=[-1, 0, 1], extend='min', aspect=40, shrink=.8)
    plt.savefig('/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/codes/Appendix/A1_corr.png')
    
def detrending(df):
    temp = df.copy()
    for var in temp.columns.values:
        print(var)
        if var == 'trend':
            pass
        else:
            model = lm.OLS(df[var],df['trend'],missing='drop')
            reg = model.fit()
            temp[var] = temp[var] - temp['trend']*reg.params[0]
    del temp['trend']
    return temp


def load_data():
    ### Set wkdir
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/Analysis'
    os.chdir(wkdir)
    ### Read price and physical data
    prices = pd.read_stata('data/transformed_data/transformed_data_prices_v19.dta')
    physicals = pd.read_stata('data/transformed_data/transformed_data_physical_v19.dta')
    sdf = pd.read_csv('data/SDF/SDF-full-sample-2020-04-13.csv').rename(columns={'Unnamed: 0':'date'})
    sdf['date'] = sdf['date'].apply(lambda x: pd.Timestamp(x))
    
    date_cols_price = [x for x in list(prices.columns.values) if ('date' in x) or ('trend' in x) or ('StikIdx' in x) or ('sdf' in x)]
    prices = prices.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(prices.columns.values) if x not in date_cols_price})
    prices = prices.rename(columns={'date_Fri':'date'})
    date_cols_physical = [x for x in list(physicals.columns.values) if ('date' in x) or ('trend' in x) or ('StikIdx' in x) or ('sdf' in x)]
    physicals = physicals.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(physicals.columns.values) if x not in date_cols_physical})
    physicals = physicals.rename(columns={'date_Tue':'date'})
    # sdf tuesday & thursday
    tues = [x for x in sdf['date'] if (int(re.findall(r'[0-9]+',str(x-sdf['date'][0]))[0])%7==4)]
    sdf_tues = sdf[sdf['date'].isin(tues)]
    sdf_tues['date'] = sdf_tues['date'].apply(lambda x: x + pd.Timedelta('3days'))
    sdf_tues.rename(columns={'E_wti_fut_tr':'sdf_fullSample'}, inplace=True)
    
    thurs = [x for x in sdf['date'] if (int(re.findall(r'[0-9]+',str(x-sdf['date'][0]))[0])%7==6)]
    sdf_thurs = sdf[sdf['date'].isin(thurs)]
    sdf_thurs['date'] = sdf_thurs['date'].apply(lambda x: x + pd.Timedelta('1days'))
    sdf_thurs.rename(columns={'E_wti_fut_tr':'sdf_fullSample'}, inplace=True)
    
    # VOL premium
    prices = pd.merge(prices, sdf_thurs, on='date', how='left')
    physicals = pd.merge(physicals, sdf_tues, on='date', how='left')
    
    # Selected Variables
    price_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                    'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI_8wk', 'trend', 'artcount', 
                    'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv', 'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp', 
                    'VIX', 'ovx_diff', 'vix_diff', 'sdf_fullSample', 'PCAsent', 'PCAfreq', 'PCAall', 'BEME', 'Mom', 'BasMom', 
                    'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']
    physical_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                    'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI_8wk', 'trend', 'artcount', 
                    'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv', 'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp', 
                    'VIX', 'ovx_diff', 'vix_diff', 'sdf_fullSample', 'PCAsent', 'PCAfreq', 'PCAall', 'BEME', 'Mom', 'BasMom', 
                    'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']
    prices = prices.loc[:,price_vars]
    physicals = physicals.loc[:,physical_vars]
    
    # Detrend
    prices = detrending(prices)
#    physicals = detrending(physicals)
    return prices

def clustering(df):
    cluster_th = 1
    
    X = df.corr().values
    d = sch.distance.pdist(X)
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    
    columns = [df.columns.tolist()[i] for i in list(np.argsort(ind))]
    df = df.reindex(columns, axis=1)
    
    unique, counts = np.unique(ind, return_counts=True)
    counts = dict(zip(unique, counts))
    
    i = 0
    j = 0
    columns = []
    for cluster_l1 in set(sorted(ind)):
        j += counts[cluster_l1]
        sub = df[df.columns.values[i:j]]
        if counts[cluster_l1]>cluster_th:        
            X = sub.corr().values
            d = sch.distance.pdist(X)
            L = sch.linkage(d, method='complete')
            ind = sch.fcluster(L, 0.2*d.max(), 'distance')
            col = [sub.columns.tolist()[i] for i in list((np.argsort(ind)))]
            sub = sub.reindex(col, axis=1)
        cols = sub.columns.tolist()
        columns.extend(cols)
        i = j
    df = df.reindex(columns, axis=1)
    return df

def main():
    prices = load_data()
    prices_clustered = clustering(prices)
#    physicals_clustered = clustering(physicals)
#    plot_corr(prices, size=12)
    plot_corr(prices_clustered, size=12)
#    plot_corr(physicals, size=12)
#    plot_corr(physicals_clustered, size=12)

if __name__=='__main__':
    main()
