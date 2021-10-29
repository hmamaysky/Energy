#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:03:22 2020

This code pick the best one and one model for each dependent variable based on RMSE.

The RMSE calculation is conducted using the other code, this is just the summary part.

This model fixed the 1-1 pair between text and baseline variables at first,
 and update the coefficients using 5-yr backward looking window as well as making prediction weekly.
 The results are rather promising because we are optimizing the RMSE using all the data.

@author: billwu
"""
import pandas as pd
import numpy as np
import pickle
import os
import maplotlib.pyplot as plt

# %%
def get_date_index(date_dict={0:('2003-04-25','2007-11-30'),
                              1:('2007-12-01','2009-06-30'),
                              2:('2009-07-01','2015-10-31'),
                              3:('2015-11-01','2020-01-31')}):
    '''
    Get date index for start and end date in each period.
    ovx_diff and sdf series starts later than other series. 
    Make sure to take care of that.
    @date_dict: the separated period we would like to test. We want to see how the models behave in those intervals
    '''
    index_dict={'all':[],'ovx':[],'sdf':[]}
    all_dates = {'all':[str(pd.Timestamp('2003-04-25')+pd.Timedelta('{} days'.format(7*x)))[:10] for x in range(1000)],
                 'sdf':[str(pd.Timestamp('2002-04-08')+pd.Timedelta('{} days'.format(7*x+7*5*52)))[:10] for x in range(1000)],
                 'ovx':[str(pd.Timestamp('2007-05-11')+pd.Timedelta('{} days'.format(7*x+7*5*52)))[:10] for x in range(1000)]}
    index_dict={'all':[], 'sdf':[], 'ovx':[]}  
    for spec in index_dict.keys():
        for dates in date_dict.values():
            if dates[0]>=min(all_dates[spec]):
                index_dict[spec].append((len([x for x in all_dates[spec] if x <= dates[0]])-1, 
                                         len([x for x in all_dates[spec] if x <= dates[1]])-1))
            else:
                index_dict[spec].append((0,0))
    return index_dict

def data_concat():
    data = pd.read_excel(files[-1])

    for i in range(len(files)-1):
        temp = pd.read_excel(files[i])
        data = pd.concat([data, temp.iloc[1:,:]])
    data.rename(columns={'Unnamed: 0': 'model'}, inplace = True)
    data=data[data['model']!='Constant']
    return data

def var_key_map(var):
    if 'sdf' in var:
        return 'sdf'
    elif 'ovx' in var:
        return 'ovx'
    else:
        return 'all'
        
def less_than_const(date_dict={0:('2003-04-25','2007-11-30'),
                              1:('2007-12-01','2009-06-30'),
                              2:('2009-07-01','2015-10-31'),
                              3:('2015-11-01','2020-01-31')}):
    '''
    This function output the models that beats the constant model and their MSE ratios against the constant.
    There are 8 dependent vars, so 8 .csv files are produced.
    For ovx, and sdf series, the time horizon is different from the regular variables, so the calculation,
    and the comparison should be separately treated.
    '''
    ### Start the main process
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/new_variables/fixed_model/oneandone/weekly/'
    wkdir += 'rmse'
    os.chdir(wkdir)
    
    ### Start Processing, get the rmse on start and end dates in the date_dict
    # all the vars
    var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot',
                'DOilVol', 'DInv', 'DProd']
    result_dict = {x:{} for x in var_list} # list to hold the results
    # Get index of the start and end dates (see comments in the get_date_index function)
    index_dict = get_date_index(date_dict)
    for rmse_p in os.listdir():
        temp_p = pickle.load(open(rmse_p,'rb'))
        print(rmse_p)
        for var in var_list:
            # for each specification in the dataframe (for each var, and if it's constant or not)
            for key in [x for x in temp_p.keys() if var in x]:
                if 'cons' in key:
                    temp_mse_list=[]
                    for date_ind in index_dict[var_key_map(rmse_p)]: 
                        if str(date_ind)=='(0,0)':
                            temp_mse_list.append(np.nan)
                        elif date_ind[0]==0:
                            temp_mse_list.append(np.power(temp_p[key][date_ind[1]],2))
                        else:
                            temp_mse_list.append((np.power(temp_p[key][date_ind[1]],2)*date_ind[1]\
                            -np.power(temp_p[key][date_ind[0]],2)*date_ind[0])\
                           /(date_ind[1]-date_ind[0]-1))        
                    result_dict[var]['{}_const'.format(rmse_p.split('_')[0])] = temp_mse_list
                else:
                    for key_t in temp_p[key].keys():
                        if ('sdf' not in rmse_p) and ('sdf' in key_t):
                            continue
                        temp_mse_list=[]
                        for date_ind in index_dict[var_key_map(rmse_p)]: 
                            if str(date_ind)=='(0,0)':
                                temp_mse_list.append(np.nan)
                            elif date_ind[0]==0:
                                temp_mse_list.append(np.power(temp_p[key][key_t][date_ind[1]],2))
                            else:
                                temp_mse_list.append((np.power(temp_p[key][key_t][date_ind[1]],2)*date_ind[1]\
                                -np.power(temp_p[key][key_t][date_ind[0]],2)*date_ind[0])\
                               /(date_ind[1]-date_ind[0]-1))        
                        result_dict[var][key_t] = temp_mse_list      
    # Further process the result_dict
    result_dict = {x:pd.DataFrame(result_dict[x]) for x in var_list}
    for keys in result_dict.keys():
        result_dict[keys].index = [x[0]+' to '+x[1] for x in date_dict.values()]
    ### Analysis for each D_VAR
    print('Get Winning Models ...')
    final_result_dict = {}
    for d_var in var_list:
        print(d_var)
        final_result_dict[d_var]={}
        result_df = result_dict[d_var]
        # Separate regular, ovx and sdf variables
        ovx_col = [x for x in result_df.columns.values if 'ovx' in x]
        sdf_col = [x for x in result_df.columns.values if 'sdf' in x and 'ovx' not in x]
        reg_col = [x for x in result_df.columns.values if x not in sdf_col and x not in ovx_col]
        mse_ovx = result_df.loc[:,ovx_col]
        mse_sdf = result_df.loc[:,sdf_col]
        mse_reg = result_df.loc[:,reg_col]
        
        # Only keep one col for constant value
        const_col = [x for x in mse_ovx.columns.values if 'const' in x][0]
        mse_ovx = mse_ovx.rename(columns={const_col:'Const'})
        mse_ovx = mse_ovx.drop([x for x in mse_ovx.columns.values if 'const' in x], axis=1)
        const_col = [x for x in mse_sdf.columns.values if 'const' in x][0]
        mse_sdf = mse_sdf.rename(columns={const_col:'Const'})
        mse_sdf = mse_sdf.drop([x for x in mse_sdf.columns.values if 'const' in x], axis=1)
        const_col = [x for x in mse_reg.columns.values if 'const' in x][0]
        mse_reg = mse_reg.rename(columns={const_col:'Const'})
        mse_reg = mse_reg.drop([x for x in mse_reg.columns.values if 'const' in x], axis=1)
        
        ### Calculate MSE ratio
        mse_ovx = mse_ovx.div(mse_ovx['Const'],axis=0).drop('Const',axis=1)
        mse_sdf = mse_sdf.div(mse_sdf['Const'],axis=0).drop('Const',axis=1)
        mse_reg = mse_reg.div(mse_reg['Const'],axis=0).drop('Const',axis=1)
        
        ### Concat 
        mse_all = pd.concat([mse_ovx.T, mse_sdf.T, mse_reg.T], axis=0)
        ### Get Ratio for each
        for periods in list(mse_all.columns.values):
            temp_df = mse_all.loc[:,[periods]].reset_index().rename(columns={'index':'model',periods:'ratio'})
            temp_df = temp_df.dropna(subset=['ratio'])
            temp_df = temp_df[temp_df.ratio<1].sort_values(by=['ratio'],ascending=False)
            temp_df['set'] = temp_df.loc[:,'model'].apply(lambda x: frozenset(x.split(', ')))
            temp_df = temp_df.drop_duplicates(subset='set')
            temp_df = temp_df.drop('set',axis=1)
            final_result_dict[d_var][periods] = temp_df
        
    return final_result_dict    

def summary_stats(win_model_dict, date_dict={0:('2003-04-25','2007-11-30'),
                                             1:('2007-12-01','2009-06-30'),
                                             2:('2009-07-01','2015-10-31'),
                                             3:('2015-11-01','2020-01-31')}):
    # Analysis for each period
    period_list = [x[0]+' to '+x[1] for x in date_dict.values()]
    for period in period_list:
        no_text=20; no_base=28
        # MAKE SURE ovx and sdf are excluded properly
        if period[:10]<=str(pd.Timestamp('2002-04-08')+pd.Timedelta('{} days'.format(7*5*52)))[:10]:
            no_base-=3
        elif period[:10]<=str(pd.Timestamp('2007-05-11')+pd.Timedelta('{} days'.format(7*5*52)))[:10]:
            no_base-=1
        # Compute quantiles and stats of the MSE ratios of the winning models for each dependent variable
        summary_df = {}
        var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot',
           'DOilVol', 'DInv', 'DProd']
        # list for specific type of models
        text_vars = ['artcount', 'entropy', 'sCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
                  'sEpg', 'fCo', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
                  'sent', 'PCAsent', 'PCAfreq','PCAall']
        text_models = []
        base_models = []
        good_models = []
        # main summary process
        for var in var_list:
            temp_ = win_model_dict[var][period]
            temp = temp_['ratio']
            stats = temp.quantile([.05, .10, .25, .5]).tolist()
            stats.extend([np.mean(temp), np.std(temp), '{}/{}'.format(len(temp), 47*23)])
            summary_df[var] = stats
            text=0;base=0
            for model in temp_.model.tolist():
                matched = [x for x in text_vars if x in model]
                if len(matched) == 2:
                    text+=1
                elif len(matched) == 1:
                    text+=1
                    base+=1
                else:
                    base+=1
            temp_good = len(temp[temp<0.95])
            text_models.append('{}/{}'.format(text,int(no_text*(no_text-1)/2+no_text*no_base)))
            base_models.append('{}/{}'.format(base,int(no_base*(no_base-1)/2+no_text*no_base)))
            good_models.append('{}/{}'.format(temp_good,len(temp)))
           
        summary_df = pd.DataFrame(summary_df).transpose()
        summary_df.rename(columns={0:'p5',1:'p10',2:'p25',3:'p50',4:'mean',5:'std',6:'N(Winning)'}, inplace=True)
        summary_df['N(Wining, <0.95)'] = good_models
        summary_df['N(Text)'] = text_models
        summary_df['N(Base)'] = base_models
        summary_df.to_csv('../mse_summary/{}_stats.csv'.format(period),index=False)
        
def summary_driver(date_dict={0:('2003-04-25','2007-11-30'),
                              1:('2007-12-01','2009-06-30'),
                              2:('2009-07-01','2015-10-31'),
                              3:('2015-11-01','2020-01-31')}):
    # Let's get started
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/new_variables/fixed_model/oneandone/weekly/'
    wkdir += 'rmse'
    os.chdir(wkdir)
    data_proc = less_than_const()
    summary_stats(data_proc,date_dict)
    
# %%
if __name__ == "__main__":
    time_dict={0:('2003-04-25','2007-11-30'),
               1:('2007-12-01','2009-06-30'),
               2:('2009-07-01','2015-10-31'),
               3:('2015-11-01','2020-01-31')}
    summary_driver(date_dict=time_dict)
    
