#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 23:04:36 2020

@author: hw2676

This file counts the total number of time that a model is selected in the stability 
test thoughout the whole time span for each predicted variable.
The list of selected model is produced by select_models.py, so the process here
is quite straightforward.
"""

import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import os
import sys
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
              'bmratio', 'mom_fut1', 'basismom', 'dxy_betas', 'cpiyr_betas', 'hp', 'liquidity', 'oi',
              'trend', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
    all_vars = text_vars + base_vars
    matrix = pd.DataFrame(0, index=all_vars, columns=all_vars)
    return matrix
    

def matrix_addition(mat, model_list):
    for i in model_list:
        i=i.split(', ')
        mat.loc[i[0], i[1]]+=1
        mat.loc[i[1], i[0]]+=1
    return mat

def matrix_process_weekly(mat, model_dict):
    model_list = list(model_dict.keys())
    return matrix_addition(mat, model_list)

   
def matrix_plot(mat,lookback=3):
    fig, ax = plt.subplots(figsize=(31,32), dpi=200)
    cax = ax.matshow(mat, cmap=plt.cm.YlOrRd)
    fig.colorbar(cax, ticks=[0,1,2,3,4,5,6,7,8], aspect=40, shrink=.6)    
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    base_vars = ['FutRet', 'xomRet', 'rdsaRet', 'bpRet', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPImom', 
              'bmratio', 'mom_fut1', 'basismom', 'dxy_betas', 'cpiyr_betas', 'hp', 'liquidity', 'oi',
              'trend', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
    all_vars = text_vars + base_vars
    ax.set_xticks(np.arange(len(all_vars)))
    ax.set_xticklabels(all_vars, rotation=90)
    ax.set_yticks(np.arange(len(all_vars)))
    bottom, top = ax.get_ylim()
    ax.hlines([19.5],-0.5,46.5,color='black')
    ax.vlines([19.5],-0.5,46.5,color='black')
#    ax.set_ylim(bottom+0.5, top-0.5)
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
    rect_1 = patches.Rectangle((mat.shape[1]+0.4, 8.9),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')
    rect_2 = patches.Rectangle((mat.shape[1]+0.4, 32.4),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')    
    rect_1.set_clip_on(False)
    rect_2.set_clip_on(False)
    ax.add_patch(rect_1)
    ax.add_patch(rect_2)
    plt.savefig('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/matrix_plots/'+
                str(lookback)+'yr_models.png')
    
def main(lookback=3):
    wkdir = '/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/processed'
    os.chdir(wkdir)
    
    selected_models = pickle.load(open(str(lookback)+'yrback_selected_models_with_coefs.p','rb'))
    matrix = matrix_generate()
    
    d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    for d_var in d_vars:
        print(d_var)
        model_list = selected_models[d_var]
        for week, model_dict in model_list.items():
            matrix = matrix_process_weekly(matrix, model_dict)
            
    # matrix_plot(matrix,lookback=lookback)
    return matrix


if __name__ == "__main__":
    # Get lookback window from command line
    # window = sys.argv[1]
    window = '3'
    # Perform whole analysis
    matrix = main(lookback=int(window))
    matrix_plot(matrix,lookback=window)
    matrix.to_csv('/user/hw2676/files/Energy/outputs/wipimom_updated/new_variables/fixed_model/time_varying_model/matrix_plots/'+
                window+'yr_models_raw_matrix.csv')