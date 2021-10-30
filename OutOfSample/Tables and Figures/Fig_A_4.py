#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 17:10:32 2020

@author: billwu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import os
import warnings
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['axes.labelweight'] = 'bold'
warnings.filterwarnings('ignore')

def matrix_generate():
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    base_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                    'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI_8wk', 'trend', 
                    'ovx_diff', 'vix_diff', 'RPsdf_rolling', 'RPsdf_growing', 'BEME', 'Mom', 'BasMom', 
                    'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']
    all_vars = text_vars + base_vars
    matrix = pd.DataFrame(0, index=all_vars, columns=all_vars)
    return matrix
    

def matrix_process(mat, series):
    temp = series.apply(lambda x: x.split(', '))
    for i in temp:
        mat.loc[i[0], i[1]]+=1
        mat.loc[i[1], i[0]]+=1
    return mat

def matrix_plot(mat):
    fig, ax = plt.subplots(figsize=(17,18), dpi=200)
    cax = ax.matshow(mat, cmap=plt.cm.YlOrRd)
    fig.colorbar(cax, ticks=[0,1,2,3,4,5,6,7,8], aspect=40, shrink=.7)    
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    base_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                 'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI', 'trend', 
                 'ovx_diff', 'vix_diff', 'sdf_rolling', 'sdf_growing', 'BE/ME', 'Mom', 'BasMom', 
                 'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']
    all_vars = text_vars + base_vars
    ax.set_xticks(np.arange(len(all_vars)))
    ax.set_xticklabels(all_vars, rotation=90)
    ax.set_yticks(np.arange(len(all_vars)))
    bottom, top = ax.get_ylim()
    ax.hlines([19.5],-0.5,47.5,color='black')
    ax.vlines([19.5],-0.5,47.5,color='black')
    ax.set_ylim(bottom, top)
    ax.set_yticklabels(all_vars)
    for i in range(mat.shape[0]):
        ax.text(mat.shape[1], i, sum(mat.iloc[i,:]), va='center', ha='center')
        for j in range(mat.shape[1]):
            c = mat.iloc[i,j]
            ax.text(j, i, str(c), va='center', ha='center')
            if i ==1:
                ax.text(j, mat.shape[0], sum(mat.iloc[:,j]), va='center', ha='center')
    ax.text(mat.shape[1]+1, 9.5, mat.iloc[:20,:].sum().sum(), va='center', ha='center')
    ax.text(mat.shape[1]+1, 33, mat.iloc[20:,:].sum().sum(), va='center', ha='center')
#    rect_1 = patches.Rectangle((mat.shape[1]+0.4, 8.9),1.2,1.2,
#                             linewidth=1,edgecolor='black',facecolor='none')
#    rect_2 = patches.Rectangle((mat.shape[1]+0.4, 32.4),1.2,1.2,
#                             linewidth=1,edgecolor='black',facecolor='none')    
#    rect_1.set_clip_on(False)
#    rect_2.set_clip_on(False)
#    ax.add_patch(rect_1)
#    ax.add_patch(rect_2)
    plt.tight_layout()
    plt.savefig('../matrix_plots/matrix_plot.png')
    
def main():
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/new_variables/fixed_model/oneandone/weekly'
    wkdir += '/results'
    os.chdir(wkdir)
    
    matrix = matrix_generate()
    files = os.listdir()
    try:
        files.remove('.DS_Store')
    except:
        pass
    for i in files:
        series = pd.read_csv(i)['model']
        matrix = matrix_process(matrix, series)
    
    matrix_plot(matrix)
    return matrix


if __name__ == "__main__":
    matrix = main()

    
