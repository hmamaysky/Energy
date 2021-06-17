#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 17:10:32 2020

@author: Hongyu Wu

This code makes the matrix plot in the paper.
The plot has all the variables on the top and on the left hand side edge and is symetric.
The text vars precede the baseline vars in the var names.
The blocks indicate the number of times that the corresponding fixed model beats the constant 
model in all 8 cases (there are 8 dependent variables thus 8 cases).
Total selection time of a variable is calculated and shown on the right edge.
There are two blocks indicating the total selection time of all the text vars and the baseline vars respectively
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
    base_vars = ['FutRet', 'xomRet', 'rdsaRet', 'bpRet', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPImom_8wk', 
              'trend', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
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
    fig, ax = plt.subplots(figsize=(15,16), dpi=200)
    cax = ax.matshow(mat, cmap=plt.cm.YlOrRd)
    fig.colorbar(cax, ticks=[0,1,2,3,4,5,6,7,8], aspect=40, shrink=.8)    
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    base_vars = ['FutRet', 'xomRet', 'rdsaRet', 'bpRet', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPImom', 
              'trend', 'vix_spx', 'ovx_cl1', 'sdf_growing', 'sdf_rolling']
    all_vars = text_vars + base_vars
    ax.set_xticks(np.arange(len(all_vars)))
    ax.set_xticklabels(all_vars, rotation=90)
    ax.set_yticks(np.arange(len(all_vars)))
    bottom, top = ax.get_ylim()
    ax.hlines([19.5],-0.5,38.5,color='black')
    ax.vlines([19.5],-0.5,38.5,color='black')
    ax.set_ylim(bottom+0.5, top-0.5)
    ax.set_yticklabels(all_vars)
    for i in range(mat.shape[0]):
        ax.text(mat.shape[1], i, sum(mat.iloc[i,:]), va='center', ha='center')
        for j in range(mat.shape[1]):
            c = mat.iloc[i,j]
            ax.text(j, i, str(c), va='center', ha='center')
            if i ==1:
                ax.text(j, mat.shape[0], sum(mat.iloc[:,j]), va='center', ha='center')
    ax.text(mat.shape[1]+1, 9.5, mat.iloc[:20,:].sum().sum(), va='center', ha='center')
    ax.text(mat.shape[1]+1, 29, mat.iloc[20:,:].sum().sum(), va='center', ha='center')
    rect_1 = patches.Rectangle((mat.shape[1]+0.4, 8.9),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')
    rect_2 = patches.Rectangle((mat.shape[1]+0.4, 28.4),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')    
    rect_1.set_clip_on(False)
    rect_2.set_clip_on(False)
    ax.add_patch(rect_1)
    ax.add_patch(rect_2)
    plt.savefig('../matrix_plots/matrix_plot.png')
    
def main():
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/Analysis/Prediction Power of textual'
    wkdir += '/outcome/model selection results/20201116 WIPImom Updates/wipimom_updated/final_codes_test/fixed_model/oneandone/weekly'
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

    
