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

# %% 0. Importing the Packages
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

# %% 1. Defining Functions

## 1.1 Function for the matrix generation
def matrix_generate():
    
    ''' This function generates the empty matrix with all the vars in column and index. '''
    
    # text vars
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    
    # baseline vars
    base_vars = ['FutRet', 'xomRet', 'rdsaRet', 'bpRet', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPIyoy', 
              'trend', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
    
    # all vars
    all_vars = text_vars + base_vars
    
    # generating the matrix
    matrix = pd.DataFrame(0, index=all_vars, columns=all_vars)
    
    # return the matrix
    return matrix
    
## 1.2 Function for matrix updating each round
def matrix_process(mat, series):
    '''
    This function update the matrix to be plotted by adding one to the 
    corresponding block if that model beats the constant model in the simulation
    for a certain LHS variable.
    
    Inputs:
        1. The matrix to be updated
        2. The series with all the winning fixed model specifications
    Output:
        1. The updated matrix
    '''
    # get all the winning models
    temp = series.apply(lambda x: x.split(', '))
    # add one to the matrix per winning model
    for i in temp:
        mat.loc[i[0], i[1]]+=1
        mat.loc[i[1], i[0]]+=1
    return mat

## 1.3 Function for Matrix plotting 
def matrix_plot(mat):
    
    ''' 
    This function plots the matrix after all the processing.
    Please specify the proper directory to save the plots 
    at the end of this function.
    '''
    
    ## 1. Draw the canvas and plot the matrix
    fig, ax = plt.subplots(figsize=(15,16), dpi=200)
    cax = ax.matshow(mat, cmap=plt.cm.YlOrRd)
    ## 2. Set the color bar
    fig.colorbar(cax, ticks=[0,1,2,3,4,5,6,7,8], aspect=40, shrink=.8)   
    ## 3. Set x and y tickers
    # x and y ticker names
    text_vars = ['artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
              'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
              'PCAsent', 'PCAfreq', 'PCAall']
    base_vars = ['FutRet', 'xomRet', 'rdsaRet', 'bpRet', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot',
              'tnote_10y', 'DFX', 'sp500Ret', 'basis', 'WIPIyoy', 
              'trend', 'vix_spx', 'ovx_cl1', 'RPsdf_growing', 'RPsdf_rolling']
    all_vars = text_vars + base_vars
    # set ticks and ticklabels
    ax.set_xticks(np.arange(len(all_vars)))
    ax.set_xticklabels(all_vars, rotation=90)
    ax.set_yticks(np.arange(len(all_vars)))
    ax.set_yticklabels(all_vars)
    ## 4. Draw a line to seperate text and baseline vars
    ax.hlines([19.5],-0.5,38.5,color='black')
    ax.vlines([19.5],-0.5,38.5,color='black')
    ## 5. Adjust the y limit toget proper plot
    # if not, the top row will only display half of the block
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom+0.5, top-0.5)
    ## 6. Annotation inside the matrix plot
    for i in range(mat.shape[0]):
        # annotation in the center of the blocks the selection time
        ax.text(mat.shape[1], i, sum(mat.iloc[i,:]), va='center', ha='center')
        for j in range(mat.shape[1]):
            c = mat.iloc[i,j]
            ax.text(j, i, str(c), va='center', ha='center')
            if i == 1:
                # annotate the sum of selected time for a var on the right hand side edge
                ax.text(j, mat.shape[0], sum(mat.iloc[:,j]), va='center', ha='center')
    # annotate the total selected time of the text and baseline vars and surrond them by a square
    ax.text(mat.shape[1]+1, 9.5, mat.iloc[:20,:].sum().sum(), va='center', ha='center')
    ax.text(mat.shape[1]+1, 29, mat.iloc[20:,:].sum().sum(), va='center', ha='center')
    # draw the square mentioned just now and add them to the canvas at the desired position
    rect_1 = patches.Rectangle((mat.shape[1]+0.4, 8.9),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')
    rect_2 = patches.Rectangle((mat.shape[1]+0.4, 28.4),1.2,1.2,
                             linewidth=1,edgecolor='black',facecolor='none')    
    rect_1.set_clip_on(False)
    rect_2.set_clip_on(False)
    ax.add_patch(rect_1)
    ax.add_patch(rect_2)
    plt.savefig('Please specify the proper directory to save the image')
    
# %% 2. Main Process    
def main():
    ## Set working directories for the files summarizing the winning fixed models (8 files, each for one LHS var)
    wkdir = 'Please Enter the Proper Directory for the results from the file "best_ratio.py"'
    os.chdir(wkdir)
    
    ## Generating the matrix
    matrix = matrix_generate()
    ## Get all the files
    files = os.listdir()
    ## Process the matrix
    for i in files:
        series = pd.read_csv(i)['model']
        matrix = matrix_process(matrix, series)
    ## Plot the matrix
    matrix_plot(matrix)
    ## Return the matrix
    return matrix

# %% 3. Executing the main process
if __name__ == "__main__":
    matrix = main()
    
    