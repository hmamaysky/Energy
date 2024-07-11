#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 22:35:02 2021

@author: billwu
"""
# %%
import pandas as pd
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt

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

def var_key_map(var):
    if 'sdf' in var:
        return 'sdf'
    elif 'ovx' in var:
        return 'ovx'
    else:
        return 'all'
    
def date_shift(date, wk=8):
    ### we forecast 8 wk ahead, so drop last 8 observations in each period to avoid overlapping
    return str(pd.Timestamp(date)+pd.Timedelta('{} days'.format(7*wk)))[:10]
        
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
    for rmse_p in [x for x in os.listdir() if x.endswith('.p')]:           
        temp_p = pickle.load(open(rmse_p,'rb'))
#        print(rmse_p)
        for var in var_list:
            # for each specification in the dataframe (for each var, and if it's constant or not)
            for key in [x for x in temp_p.keys() if var in x]:
                if 'cons' in key:
                    temp_mse_list=[]
                    ii=0
                    for date_ind in index_dict[var_key_map(rmse_p)]: 
                        ii+=1
                        if ii == len(index_dict[var_key_map(rmse_p)]):
                            temp_mse_list.append((np.power(temp_p[key][date_ind[1]],2)*date_ind[1]\
                            -np.power(temp_p[key][date_ind[0]],2)*date_ind[0])\
                           /(date_ind[1]-date_ind[0]-1))                              
                        elif str(date_ind)=='(0, 0)':
                            temp_mse_list.append(np.nan)
                        elif date_ind[0]==0:
                            temp_mse_list.append(np.power(temp_p[key][date_ind[1]-8],2))
                        else:
                            temp_mse_list.append((np.power(temp_p[key][date_ind[1]-8],2)*date_ind[1]\
                            -np.power(temp_p[key][date_ind[0]-8],2)*date_ind[0])\
                           /(date_ind[1]-date_ind[0]-1))        
                    result_dict[var]['{}_const'.format(rmse_p.split('_')[0])] = temp_mse_list
                else:
                    for key_t in temp_p[key].keys():
                        if ('sdf' not in rmse_p) and ('sdf' in key_t):
                            continue
                        if ('ovx' not in rmse_p) and ('ovx' in key_t):
                            continue
                        temp_mse_list=[]
                        ii=0
                        for date_ind in index_dict[var_key_map(rmse_p)]: 
                            ii+=1
                            if ii == len(index_dict[var_key_map(rmse_p)]):
                                temp_mse_list.append((np.power(temp_p[key][key_t][date_ind[1]],2)*date_ind[1]\
                                -np.power(temp_p[key][key_t][date_ind[0]],2)*date_ind[0])\
                               /(date_ind[1]-date_ind[0]-1))   
                            elif str(date_ind)=='(0, 0)':
                                temp_mse_list.append(np.nan)
                            elif date_ind[0]==0:
                                temp_mse_list.append(np.power(temp_p[key][key_t][date_ind[1]-8],2))
                            else:
                                temp_mse_list.append((np.power(temp_p[key][key_t][date_ind[1]-8],2)*date_ind[1]\
                                -np.power(temp_p[key][key_t][date_ind[0]-8],2)*date_ind[0])\
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
            stats = ['{:.4f}'.format(x) for x in stats]
            stats.extend(['{:.4f}'.format(np.mean(temp)), '{:.4f}'.format(np.std(temp)), '{}/{}'.format(len(temp), int((no_text+no_base)*(no_text+no_base-1)/2))])
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
        summary_df['N(Nontext)'] = base_models
        summary_df.to_csv('../mse_summary/{}_stats.csv'.format(period))
       
def transition_matrix(win_model_dict, date_dict={0:('2003-04-25','2007-11-30'),
                                                 1:('2007-12-01','2009-06-30'),
                                                 2:('2009-07-01','2015-10-31'),
                                                 3:('2015-11-01','2020-01-31')},
                                      mode='all'):
    def positive_models(date0,mode='all'):
        '''Helper function to get the set of candidate models'''
        all_dates = {'all':'2003-04-25',
                     'sdf':str(pd.Timestamp('2002-04-08')+pd.Timedelta('{} days'.format(7*5*52)))[:10],
                     'ovx':str(pd.Timestamp('2007-05-11')+pd.Timedelta('{} days'.format(7*5*52)))[:10]}

        all_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                    'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI_8wk', 'trend', 'artcount', 
                    'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv', 'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp', 
                    'VIX', 'ovx_diff', 'vix_diff', 'RPsdf_rolling', 'RPsdf_growing', 'PCAsent', 'PCAfreq', 'PCAall', 'BEME', 'Mom', 'BasMom', 
                    'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt', 'sent']
        base_vars = ['FutRet', 'DSpot', 'DOilVol', 'xomRet', 'bpRet', 'rdsaRet',
                    'OilVol', 'DInv', 'DProd', 'tnote_10y', 'DFX', 'sp500Ret', 'StkIdx', 'basis', 'WIPI_8wk', 'trend', 
                    'VIX', 'ovx_diff', 'vix_diff', 'RPsdf_rolling', 'RPsdf_growing', 'BEME', 'Mom', 'BasMom', 
                    'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']
        if date0>=all_dates['ovx']:
            pass
        elif date0>=all_dates['sdf']:
            all_vars.remove('ovx_diff')
            base_vars.remove('ovx_diff')
        else:
            all_vars.remove('RPsdf_rolling')
            all_vars.remove('RPsdf_growing')
            base_vars.remove('RPsdf_rolling')
            base_vars.remove('RPsdf_growing')
        if mode =='all':                
            set1 = set('{}, {}'.format(x,y) for x in all_vars for y in all_vars if x!=y)
        elif mode =='base':
            set1 = set('{}, {}'.format(x,y) for x in base_vars for y in base_vars if x!=y)
        else:
            set1 = set('{}, {}'.format(x,y) for x in all_vars for y in all_vars if x!=y)
            set2 = set('{}, {}'.format(x,y) for x in base_vars for y in base_vars if x!=y)
            set1 = set1-set2
        return  set(frozenset(x.split(', ')) for x in set1)
    text_vars= ['artcount', 'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv', 'sEpg', 'fEpg', 'sBbl', 'fBbl', 
                'sRpc', 'fRpc', 'sEp', 'fEp', 'PCAsent', 'PCAfreq', 'PCAall']        
    ###### Start the main process
    ### All Models ###
    result_dict={}
    pospos_dict={}
    for d_vars in win_model_dict.keys():
#        print(d_vars)
        result_dict[d_vars] = {}
        pospos_dict[d_vars] = {}
        for i in range(len(date_dict)-1):
            result_dict[d_vars]['({}) to ({})'.format(i+1,i+2)]=np.zeros((2,2))
            if mode =='all':
                prev_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i][0], date_dict[i][1])]['model']))
                new_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i+1][0], date_dict[i+1][1])]['model']))
            elif mode =='text':
                prev_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i][0], date_dict[i][1])]['model']) if max([y in x for y in text_vars])==1)
                new_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i+1][0], date_dict[i+1][1])]['model']) if max([y in x for y in text_vars])==1)
            elif mode =='base':
                prev_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i][0], date_dict[i][1])]['model']) if max([y in x for y in text_vars])==0)
                new_models = set(frozenset(x.split(', ')) for x in list(win_model_dict[d_vars]['{} to {}'.format(date_dict[i+1][0], date_dict[i+1][1])]['model']) if max([y in x for y in text_vars])==0)          
            candidate_models = positive_models(date_dict[i][0],mode=mode)
            pos_pos = [x for x in prev_models if x in new_models]
            pos_neg = [x for x in prev_models if x not in new_models]
            neg_set = [x for x in candidate_models if x not in prev_models]
            neg_pos = [x for x in neg_set if x in new_models]
            neg_neg = [x for x in neg_set if x not in new_models]
            result_dict[d_vars]['({}) to ({})'.format(i+1,i+2)][0,0] = len(pos_pos)/len(prev_models)
            result_dict[d_vars]['({}) to ({})'.format(i+1,i+2)][0,1] = len(pos_neg)/len(prev_models)
            result_dict[d_vars]['({}) to ({})'.format(i+1,i+2)][1,0] = len(neg_pos)/len(neg_set)
            result_dict[d_vars]['({}) to ({})'.format(i+1,i+2)][1,1] = len(neg_neg)/len(neg_set)
            pospos_dict[d_vars]['({}) to ({})'.format(i+1,i+2)] = [str(x) for x in pos_pos]

    ###### Get models that constantly work
    always_win = {}
    for d_vars in pospos_dict.keys():
        set_temp=set()
        i=0
        for model_set in pospos_dict[d_vars].values():
            if i==0:
                set_temp = set(model_set)
                i+=1
            else:
                set_temp = set_temp.intersection(model_set)
        always_win[d_vars] = list(set_temp)
    return result_dict, (pospos_dict, always_win) 

def plot_transmatrix(trans_mat_all, trans_mat_text, trans_mat_base, date_dict={0:('2003-04-25','2007-11-30'),
                                          1:('2007-12-01','2009-06-30'),
                                          2:('2009-07-01','2015-10-31'),
                                          3:('2015-11-01','2020-01-31')}):
    plt.ioff()
     
    for d_vars in trans_mat_all.keys():
        fig,ax = plt.subplots(2, 2,figsize=(12,6),dpi=200) 
        # Plots
        ax[0,0].plot([x for x in trans_mat_all[d_vars].keys()], [x[1,1] for x in trans_mat_all[d_vars].values()],linewidth=2,label='All')
        ax[0,0].plot([x for x in trans_mat_text[d_vars].keys()], [x[1,1] for x in trans_mat_text[d_vars].values()],'--',linewidth=2,label='Text')
        ax[0,0].plot([x for x in trans_mat_base[d_vars].keys()], [x[1,1] for x in trans_mat_base[d_vars].values()],'o--',linewidth=2,label='Nontext')
        ax[0,0].set_title('Negative to Negative',fontsize=12)
        ax[0,0].legend(fontsize=10)
        ax[0,0].set_ylim(0,1)       
#        ax[0,0].set_xlabel('; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][:3])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][3:6])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][6:]))
        
        ax[0,1].plot([x for x in trans_mat_all[d_vars].keys()], [x[1,0] for x in trans_mat_all[d_vars].values()],linewidth=2,label='All')
        ax[0,1].plot([x for x in trans_mat_text[d_vars].keys()], [x[1,0] for x in trans_mat_text[d_vars].values()],'--',linewidth=2,label='Text')
        ax[0,1].plot([x for x in trans_mat_base[d_vars].keys()], [x[1,0] for x in trans_mat_base[d_vars].values()],'o--',linewidth=2,label='Nontext')
        ax[0,1].set_title('Negative to Positive',fontsize=12)
        ax[0,1].legend(fontsize=10)
        ax[0,1].set_ylim(0,1)   
#        ax[0,1].set_xlabel('; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][:3])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][3:6])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][6:]))
        
        ax[1,0].plot([x for x in trans_mat_all[d_vars].keys()], [x[0,1] for x in trans_mat_all[d_vars].values()],linewidth=2,label='All')
        ax[1,0].plot([x for x in trans_mat_text[d_vars].keys()], [x[0,1] for x in trans_mat_text[d_vars].values()],'--',linewidth=2,label='Text')
        ax[1,0].plot([x for x in trans_mat_base[d_vars].keys()], [x[0,1] for x in trans_mat_base[d_vars].values()],'o--',linewidth=2,label='Nontext')
        ax[1,0].set_title('Positive to Negative',fontsize=12)
        ax[1,0].legend(fontsize=10)
        ax[1,0].set_ylim(0,1)   
#        ax[1,0].set_xlabel('; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][:3])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][3:6])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][6:]))
        
        ax[1,1].plot([x for x in trans_mat_all[d_vars].keys()], [x[0,0] for x in trans_mat_all[d_vars].values()],linewidth=2,label='All')
        ax[1,1].plot([x for x in trans_mat_text[d_vars].keys()], [x[0,0] for x in trans_mat_text[d_vars].values()],'--',linewidth=2,label='Text')
        ax[1,1].plot([x for x in trans_mat_base[d_vars].keys()], [x[0,0] for x in trans_mat_base[d_vars].values()],'o--',linewidth=2,label='Nontext')
        ax[1,1].set_title('Positive to Positive',fontsize=12)
        ax[1,1].legend(fontsize=10)
        ax[1,1].set_ylim(0,1)   
#        ax[1,1].set_xlabel('; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][:3])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][3:6])+
#                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][6:]))
        
        # Format
        fig.suptitle('{}\nTransition Probabilities of Fixed Model\nPositive: Winning Models; Negative: Beaten by Constant Model'.format(d_vars),fontsize=12)
        fig.supxlabel('; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][:3])+
                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][3:6])+
                      '\n'+'; '.join(['({}) {} to {}'.format(x+1,date_dict[x][0],date_dict[x][1]) for x in date_dict.keys()][6:]),fontsize=10)
    
        fig.tight_layout()
        fig.subplots_adjust(top=0.84)
        plt.savefig('../mse_summary/{}_trans_mat.png'.format(d_vars))

def summary_detail(win_model_dict, date_dict={0:('2003-04-25','2007-11-30'),
                                             1:('2007-12-01','2009-06-30'),
                                             2:('2009-07-01','2015-10-31'),
                                             3:('2015-11-01','2020-01-31')}):         
    # helper functions
    def missing_ovx_sdf(df):
        '''Apply this to countings in periods, assign -1 if data is not available yet'''
        all_dates = {'all':str(pd.Timestamp('2003-04-25'))[:10],
                     'sdf':str(pd.Timestamp('2002-04-08')+pd.Timedelta('{} days'.format(7*5*52)))[:10],
                     'ovx':str(pd.Timestamp('2007-05-11')+pd.Timedelta('{} days'.format(7*5*52)))[:10]}
        sdf_start = [n for n,i in enumerate(date_dict.values()) if i[1]>all_dates['sdf']][1]+1
        ovx_start = [n for n,i in enumerate(date_dict.values()) if i[1]>all_dates['ovx']][1]+1  
        if 'sdf' in str(df['set']):
            for j in range(sdf_start-1):
                df[str(j+1)]=-1
        elif 'ovx' in str(df['set']):
            for j in range(ovx_start-1):
                df[str(j+1)]=-1    
        return df
        
    # Analysis for each period
    text_vars = ['artcount', 'entropy', 'sCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
                  'sEpg', 'fCo', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp',
                  'sent', 'PCAsent', 'PCAfreq','PCAall']
    var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot',
           'DOilVol', 'DInv', 'DProd']
    period_list = [x[0]+' to '+x[1] for x in date_dict.values()]
    result_dict={}
    for var in var_list:
        dict_temp = win_model_dict[var]
        result_df = pd.DataFrame()
        i=1
        for key, value in dict_temp.items():
            value[str(i)]=1
            result_df = pd.concat([result_df, value], axis=0)
            i+=1
        result_df = result_df.fillna(0)
        
        ### Important step: combine {A, B} and {B, A}
        result_df['set'] = result_df.loc[:,'model'].apply(lambda x: frozenset(x.split(', ')))
        period_counts = result_df.groupby('set')['1','2','3','4','5','6','7'].sum().applymap(lambda x:x/abs(x) if x<0 else x).reset_index()     
        period_counts = period_counts.apply(lambda x:missing_ovx_sdf(x), axis=1)
        result_df = result_df.drop(['1','2','3','4','5','6','7'],axis=1)
        result_df = pd.merge(result_df,period_counts, on='set') ## metge the sum back
#        result_df['text'] = result_df.apply(lambda x: 1 if max([t in x['model'] for t in text_vars])==1 else 0,axis=1)
        result_temp= consecutive_months(result_df)
        result_df['win_max_consec'] = result_temp['win_max_consec']
        result_df['win_period'] = result_temp['win_period']
#        result_df = result_df.drop(['1','2','3','4','5','6','7','ratio'],axis=1)
        rename_dict = {str(x+1):str(y) for x,y in date_dict.items()}
        result_df = result_df.rename(columns=rename_dict)
        result_df = result_df.drop_duplicates(subset='set')
        result_df = result_df.drop('set',axis=1)
        result_df['var_1'] = result_df['model'].apply(lambda x:x.split(', ')[0])
        result_df['var_1_text'] = result_df['var_1'].apply(lambda x:1 if x in text_vars else 0)
        result_df['var_2'] = result_df['model'].apply(lambda x:x.split(', ')[1])
        result_df['var_2_text'] = result_df['var_2'].apply(lambda x:1 if x in text_vars else 0)
        result_df = result_df.drop(['ratio','model'],axis=1)
        result_dict[var] = result_df
        result_df.to_csv('/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/new_variables/fixed_model/oneandone/weekly/subperiod_details/{}.csv'.format(var)
                        ,index=False)
        
    return result_dict

def consecutive_months(df_temp):
    '''
    Find consecutive months for results in summary_detail
    '''
    df = df_temp.copy()
    df = df.drop(['ratio','model', 'set'],axis=1).T
    df = df.replace(-1,0)
    
    a = df != 0
    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df2 = df1.max()
    periods = []
    for i in range(df1.shape[1]):
        df_temp = df1.iloc[:,i]
        df_temp = df_temp[df_temp == df2[i]]
        max_temp = df2[i]
        df_index = list(df_temp.index)
        if 'set' in df_index:
            df_index.remove('set')
        periods.append('; '.join(['{}->{}'.format(int(x)-int(max_temp)+1,int(x)) for x in df_index]))
    
    df_t = df_temp.copy()
    df_t['win_max_consec'] = df2
    df_t['win_period'] = pd.Series(periods)
    
    return df_t

        
def summary_driver(date_dict={0:('2003-04-25','2007-11-30'),
                              1:('2007-12-01','2009-06-30'),
                              2:('2009-07-01','2015-10-31'),
                              3:('2015-11-01','2020-01-31')}):
    # Let's get started
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/Energy Project/revision/new_variables/fixed_model/oneandone/weekly/'
    wkdir += 'rmse'
    os.chdir(wkdir)
    data_proc = less_than_const(date_dict=date_dict)
    details = summary_detail(data_proc,date_dict=date_dict)
    summary_stats(data_proc,date_dict=date_dict)
    results_dict={}
    transmatrix_ls = []
    for mode in ['all', 'base', 'text']:
        print(mode)
        transition_mat, pospos = transition_matrix(data_proc,date_dict=date_dict,mode=mode)   
        transmatrix_ls.append(transition_mat)
        results_dict[mode+'_pospos'] = pospos[0]
        results_dict[mode+'_alwayswin'] = pospos[1]
    plot_transmatrix(transmatrix_ls[0],transmatrix_ls[1],transmatrix_ls[2],date_dict=date_dict)
    return details, data_proc
# %%
if __name__ == "__main__":
    # These time intervals are for subsample analysis, for full sample (Table A IX), just pass in the whole sample period as one interval.
    time_dict={0:('2003-04-25','2005-03-31'),
               1:('2005-04-01','2007-11-30'),
               2:('2007-12-01','2009-06-30'),
               3:('2009-07-01','2012-02-29'),
               4:('2012-03-01','2014-10-31'),
               5:('2014-11-01','2017-06-30'),
               6:('2017-07-01','2020-01-31')}
    summary_details, winning_models = summary_driver(date_dict=time_dict)
