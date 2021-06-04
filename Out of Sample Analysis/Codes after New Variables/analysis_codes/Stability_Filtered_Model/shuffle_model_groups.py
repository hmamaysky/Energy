#!/apps/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 21:10:13 2020

@author: Hongyu Wu

Description:
------------
    This code shuffles and separates the entire fixed models set into 30 lists. 
    The lists will later be passed to the 'generate_model_coefs.py' file
    in parallel on the grid to speed up the whole calculation process.
Note:
-----
    Please modify the last line in the main function to save the result in proper directory.
"""

# %% Importing packages
from itertools import combinations
import random
import pickle

# %% Defining Functions
def all_fixed_model(var_list, r=2):
    '''
    Inputs:
        1. var_list: list of all the variables
        2. # of RHS vars in the model, set default to 2
    Output:
        1. a dictionary with 30 keys referring to a separate part of the list of all the models
    '''
    # Get list of all the combination of r
    all_list = list(combinations(var_list, r))
    # Shuffle the list
    random.shuffle(all_list)
    # Save the Shuffled results
    shuffled_model_list=dict()
    for i in range(30):
        shuffled_model_list[str(i+1).zfill(2)]=all_list[i::30]
    return shuffled_model_list

# %% Main Function
def main():
    # Var list all the text and baseline vars
    vars_list = ["FutRet", "xomRet", "bpRet", "rdsaRet", "DOilVol",
                 "OilVol", "DInv", "DProd", "DSpot", "tnote_10y",
                 "DFX", "sp500Ret", "basis", "WIPImom_8wk", "trend", "RPsdf_growing",
                 "RPsdf_rolling", "vix_diff", "ovx_diff",
                 'BEME', 'Mom', 'BasMom', 'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt',
                 'artcount', 'entropy', 'sent', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
                 'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
                 'PCAfreq','PCAsent', 'PCAall']
    # Generate random list by listing all the combinations of 2 and 
    # save the result for further usage in parallelizing the coefficient calculation
    shuffled_list = all_fixed_model(vars_list,2)
    # Please change directory before running the code
    pickle.dump(shuffled_list, open('/user/hw2676/code/Energy/Analysis/New_Var/1.0Version/analysis_codes/Stability_Filtered_Model/shuffled_model_list.p','wb'))

# %% Main Process
if __name__ == '__main__':
    main()
    