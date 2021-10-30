#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 13:31:47 2020

@author: billwu
"""

import pandas as pd
import numpy as np
import os
import warnings
import re
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['axes.labelweight'] = 'bold'
warnings.filterwarnings('ignore')

wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/'
wkdir += 'new_variables/'

### Please Change wkdir, frequency, LHS and fig.suptitle Accordingly
### aspect=24 for monthly, =48 for weekly
### extent length *0.5 for monthly, *0.25 for weekly (line 52) 
### xlabel[::20] for monthly, [::40] for weekly (line 53)
### set_xticks * 0.5 for monthly, * 0.25 for weekly (line 59)
### no rotation (plt.setp) for monthly (line 61-62)
wkdir += 'parsimonious/oneandone/'
frequency = 'weekly'
wkdir += frequency
os.chdir(wkdir)
LHS ='8wk'
text_1and1 = pd.read_csv('text_vars_' + frequency + '_' + LHS +'.csv').rename(columns={'Unnamed: 0':'date'}).set_index('date').transpose().reset_index().rename(columns={'index':'date'})
base_1and1 = pd.read_csv('base_vars_' + frequency + '_' + LHS +'.csv').rename(columns={'Unnamed: 0':'date'}).set_index('date').transpose().reset_index().rename(columns={'index':'date'})
base_1and1= base_1and1.replace(['BEME','basismom','mom_fut1','dxy_betas','cpiyr_betas','hp','liquidity','oi','WIPImom_8wk','ovx_cl1','vix_spx','StikIdx'],
                               ['BE/ME','BasMom','Mom','DolBeta','InflaBeta','HedgPres','liquidity','OpenInt','WIPI','ovx_diff','vix_diff','StkIdx'])

#text_2and2 = pd.read_csv('text_vars_2and2.csv').rename(columns={'Unnamed: 0':'date'})
#base_2and2 = pd.read_csv('base_vars_2and2.csv').rename(columns={'Unnamed: 0':'date'})

d_vars = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']

### Plots for RHS selections
def plotting_base(result_set):
    oneandone_text = {}
    for i in d_vars:
        sss=pd.get_dummies(result_set[i])
        sss.index = base_1and1['date']
        oneandone_text[i]=sss

    ### First four vars
    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
    for i in np.arange(4):
        sss = oneandone_text[d_vars[i]]
        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
        xlabel = sss.index[::40]
        ylabel = [re.sub("_8wk|RP","",x) for x in sss.columns.values]    
        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
        ax[i].set_yticklabels(ylabel, fontsize=15)
        ax[i].set_title(d_vars[i], fontsize=20)
        
    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
    ax[-1].set_xticklabels(xlabel, fontsize=15)
    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")
    fig.suptitle('1-and-1 Model: Selected Nontext Variables (LHS: '+', '.join(d_vars[:4])+ ')', fontsize=30)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    plt.savefig('vars plot M/'+LHS+'1.png',bbox_inches = 'tight',
                pad_inches = 0.5)
#    plt.show()
    
    ### Next four vars
    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
    for i in np.arange(4)+4:
        sss = oneandone_text[d_vars[i]]
        i-=4
        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
        xlabel = sss.index[::40]
        ylabel = [re.sub("_8wk|RP","",x) for x in sss.columns.values]    
        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
        ax[i].set_yticklabels(ylabel, fontsize=15)
        ax[i].set_title(d_vars[i+4], fontsize=20)
        
    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
    ax[-1].set_xticklabels(xlabel, fontsize=15)
    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")   
    fig.suptitle('1-and-1 Model: Selected Nontext Variables (LHS: '+', '.join(d_vars[4:8])+ ')', fontsize=30)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    plt.savefig('vars plot M/'+LHS+'2.png',bbox_inches = 'tight',
                pad_inches = 0.5)
#    plt.show()
    
#    ### Next four vars
#    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
#    for i in np.arange(4)+8:
#        sss = oneandone_text[d_vars[i]]
#        i-=8
#        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
#        xlabel = sss.index[::40]
#        ylabel = [x for x in sss.columns.values]    
#        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
#        ax[i].set_yticklabels(ylabel, fontsize=15)
#        ax[i].set_title(d_vars[i+8], fontsize=20)
#        
#    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
#    ax[-1].set_xticklabels(xlabel, fontsize=15)
#    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
#                 rotation_mode="anchor")   
#    fig.suptitle('1-and-1 Model: Selected Baseline Variables (LHS: '+', '.join(d_vars[8:])+ ')', fontsize=30)
#    fig.tight_layout()
#    fig.subplots_adjust(top=0.93)
#    plt.savefig('vars plot M/'+LHS+'3.png',bbox_inches = 'tight',
#                pad_inches = 0.5)
#    plt.show()
    
def plotting_text(result_set):
    oneandone_text = {}
    for i in d_vars:
        sss=pd.get_dummies(result_set[i])
        sss.index = base_1and1['date']
        oneandone_text[i]=sss

    ### First four vars
    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
    for i in np.arange(4):
        sss = oneandone_text[d_vars[i]]
        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
        xlabel = sss.index[::40]
        ylabel = [re.sub("_8wk|RP","",x) for x in sss.columns.values]    
        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
        ax[i].set_yticklabels(ylabel, fontsize=15)
        ax[i].set_title(d_vars[i], fontsize=20)
        
    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
    ax[-1].set_xticklabels(xlabel, fontsize=15)
    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")
    fig.suptitle('1-and-1 Model: Selected Text Variables (LHS: '+', '.join(d_vars[:4])+ ')', fontsize=30)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    plt.savefig('vars plot M/'+LHS+'4.png',bbox_inches = 'tight',
                pad_inches = 0.5)
#    plt.show()
    
    ### Next four vars
    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
    for i in np.arange(4)+4:
        sss = oneandone_text[d_vars[i]]
        i-=4
        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
        xlabel = sss.index[::40]
        ylabel = [re.sub("_8wk|RP","",x) for x in sss.columns.values]    
        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
        ax[i].set_yticklabels(ylabel, fontsize=15)
        ax[i].set_title(d_vars[i+4], fontsize=20)
        
    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
    ax[-1].set_xticklabels(xlabel, fontsize=15)
    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")   
    fig.suptitle('1-and-1 Model: Selected Text Variables (LHS: '+', '.join(d_vars[4:8])+ ')', fontsize=30)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    plt.savefig('vars plot M/'+LHS+'5.png',bbox_inches = 'tight',
                pad_inches = 0.5)
#    plt.show()
    
#    ### Next four vars
#    fig, ax=plt.subplots(4, 1, figsize=(30,20), sharex = True, dpi=200)
#    for i in np.arange(4)+8:
#        sss = oneandone_text[d_vars[i]]
#        i-=8
#        im=ax[i].imshow(sss.transpose(),cmap='YlGnBu',extent=[0,sss.shape[0]*0.25,0,sss.shape[1]], alpha=0.8, aspect=48/sss.shape[1])
#        xlabel = sss.index[::40]
#        ylabel = [x for x in sss.columns.values]    
#        ax[i].set_yticks(np.arange(len(ylabel))+0.5)    
#        ax[i].set_yticklabels(ylabel, fontsize=15)
#        ax[i].set_title(d_vars[i+8], fontsize=20)
#        
#    ax[-1].set_xticks(np.arange(0,sss.shape[0]*0.25,10))
#    ax[-1].set_xticklabels(xlabel, fontsize=15)
#    plt.setp(ax[-1].get_xticklabels(), rotation=30, ha="right",
#                 rotation_mode="anchor")   
#    fig.suptitle('1-and-1 Selected Text Variables (LHS: '+', '.join(d_vars[8:])+ ')', fontsize=30)
#    fig.tight_layout()
#    fig.subplots_adjust(top=0.93)
#    plt.savefig('vars plot M/'+LHS+'6.png',bbox_inches = 'tight',
#                pad_inches = 0.5)
#    plt.show()
    
    
plotting_base(base_1and1)
plotting_text(text_1and1)
