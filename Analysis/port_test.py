import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
__code_dir__ = f'c:/users/{os.getenv("USERNAME")}/code/energy/'

class PortEngine:

    def __init__(self):
        self.fname_map = {'text':'0-1','nontext':'1-0','both':'1-1'}
        self.select_end = '2009-11-30'
        self.eval_start = '2010-01-01'

        ## Shell PLC trades on Euronext Amsterdam; BP is the ADR return; FutRet is return on
        ## USO ETF; XOM is from NYSE (the data are from Bloomberg). These are daily returns
        ## including dividends. Gb03 is the 3-month T-bill rate. All data are in percent.
        fname = __code_dir__+'Analysis/support/bbg_energy_with_divd_rets.csv'
        print('Reading',fname)
        self.mktd = pd.read_csv(fname,parse_dates=['Dates'],index_col=0)
        self.dpy = 365.25/self.mktd.index.to_series().diff().dt.days.mean()
        assert 258 < self.dpy < 262 ## assert returns are daily
        
        ## line up with fwd returns & double check timing
        self.rets_fwd1d = self.mktd.shift(-1)
        assert self.rets_fwd1d.xomRet['2025-05-08'] == self.mktd.xomRet['2025-05-09']
        assert self.rets_fwd1d.FutRet['2023-02-20'] == self.mktd.FutRet['2023-02-21']
        self.rets_fwd1d['gb03'] = self.mktd.gb03/self.dpy ## no shifts to the 3-month T-bill rate
        
    def read_lasso_file(self,type,var):
        '''
        Note: The first forecast in these files is from 2003-05-09 (Fri). This must be because we
        need three observations to have a tercile classification and with a 5-year lookback these
        must be the last Friday in April 2023, then Fri 2003-05-02 and 2003-05-09.
        '''        

        assert type in self.fname_map.keys()
        
        ## read the data
        fname = __code_dir__ + 'outofsample/res_Forward_lasso/' + self.fname_map[type] + '_' + var + '.csv'
        print('Reading',fname)
        return pd.read_csv(fname,parse_dates=True,index_col=0)

    def __repr__(self):

        repstr = ''
        
        for el in ['fname_map','select_end','eval_start','dpy']:
            repstr += f'{el}: {getattr(self,el)}\n'

        for el in ['mktd','rets_fwd1d']:
            repstr += f'{el}: {getattr(self,el).shape}\n'
            
        return repstr
    
    ########## portfolio testing ##########

    def port_test(self,type,signal_var,retser,thresh=0,signal_type='pred',show_plots=True,saveout=False):
        '''
        type -- which signal to use: text, nontext, both
        signal_var -- which variable's lasso model to use, e.g. FutRet, DSpot, etc.
        retser -- which market to trade, one of FutRet, xomRet, rdsaRet, and bpRet
        thresh -- threshold of forecast for taking a position
        signal_type -- one of pred, mean, true (forward looking) as the signal to use in the
                       trading simulation 
        show_plots -- do the visualization or not
        saveout -- save the output to a file?
        '''

        ## sanity checks
        assert retser in self.mktd.columns
        
        ## get the lasso signals file
        thed = self.read_lasso_file(type,signal_var)
        level = 100 if signal_var == 'FutRet' else 0 ## levels across series
        usesd = thed.true.std()
        assert (level-0.2*usesd) <= thed.true.mean() <= (level+0.2*usesd) ## make sure have the correct level
        
        ## get the excess returns and then calculate the trading weights
        idx = (self.rets_fwd1d.index >= thed.index[0]) & (self.rets_fwd1d.index <= thed.index[-1])
        portd = pd.DataFrame(self.rets_fwd1d[retser][idx] - self.rets_fwd1d.gb03[idx],columns=[retser])
        portd['weights'] = np.nan

        ## get the signal
        portd['signal'] = thed[signal_type] - level

        ## print status
        weight_AR = 0.
        print(f'Portfolio simulation for {retser} using {type} model for {signal_var} AR={weight_AR}')

        ## set the weights from the selected signals
        for ii,tt in enumerate(portd.index):

            ## if date exists in the signal data, then set the weights from it
            if not np.isnan(portd.signal[tt]):
                if portd.signal[tt] > +thresh:
                    portd.loc[tt,'weights'] = +1
                elif portd.signal[tt] < -thresh:
                    portd.loc[tt,'weights'] = -1

            ## after first period, weight persists at AR=0.9 plus add update
            if np.isnan(portd.loc[tt,'weights']):
                if ii > 0:
                    portd.loc[tt,'weights'] = portd.weights.iloc[ii-1]
                else:
                    portd.loc[tt,'weights'] = 1 ## if not signal, start at weight of 1

        ## calculate each day's portfolio return
        assert not portd.weights.isna().any()
        port_rets = portd.weights * portd[retser]

        ## get Sharpe ratio
        sr = lambda ser: ser.mean()/ser.std()*np.sqrt(self.dpy)
        print(f'Long rets SR = {sr(portd[retser])}  Port rets SR = {sr(port_rets)}')

        if show_plots:
        
            ## get the underlying return series
            jnt_rets = pd.DataFrame({retser:portd[retser],'port':port_rets})
            cum_rets = (1+jnt_rets/100).cumprod()
            cum_rets.columns = [f'{el}: {cum_rets[el].iloc[-1]:.2f}  SR: {sr(jnt_rets[el]):.3f}'
                                for el in jnt_rets.columns]

            ## show the returns
            fig, axs = plt.subplots(2,1,figsize=(8,6))
            cum_rets.plot(title=f'Cumulative excess returns for {retser} using {signal_var}\n' + \
                          f'signal "{type}" and signal type "{signal_type}"',
                          ax=axs[0],xlabel='')

            ## show the weights
            portd[['weights','signal']].ffill().clip(upper=2,lower=-2).plot(ax=axs[1],xlabel='')
            axs[1].axhline(0,color='lightgrey',linestyle='--')
            num_years = (portd.weights.index[-1] - portd.weights.index[0]).days / 365.25
            num_trades = portd.weights.diff()[portd.weights.diff() != 0].shape[0]
            axs[1].set_title(f'Number of trades per year = {num_trades/num_years:.1f}')
            plt.tight_layout()

            if saveout:
                fname = __code_dir__ + f'Analysis/results/port-test-{type}-{signal_var}-{retser}.pdf'
                print('Saving figure to',fname)
                plt.savefig(fname,bbox_inches='tight')
            
        return thed, portd.weights, sr(portd[retser]), sr(port_rets)

    def port_test_matrix(self,signal_type,saveout):
        '''
        signal_type -- the signal (mean, pred, true -- forward looking) to use in trade simulation
        saveout -- save the output?
        '''
        
        forecast_vars = ['FutRet','bpRet','rdsaRet','xomRet','DSpot']
        dep_vars = ['FutRet','bpRet','rdsaRet','xomRet']

        ## get results
        res = pd.DataFrame(np.nan,index=forecast_vars,columns=dep_vars)
        resl = pd.Series(np.nan,index=dep_vars)
        for forecast in forecast_vars:
            for dep in dep_vars:
                _, _, srl, sr = self.port_test('text',forecast,dep,thresh=0,
                                               signal_type=signal_type,show_plots=False)
                resl[dep] = srl
                res.loc[forecast,dep] = sr

        ## display results
        res.loc['Long'] = resl
        print(res)

        ## f'Signal x Return analysis using signal type {signal_type}',
        plt.figure()
        ax = sns.heatmap(res,annot=True,fmt='.3f',vmin=-0.2,vmax=0.8,cmap='RdYlGn')
        ax.set_ylabel('Signal')
        ax.set_xlabel('Return series')
        ax.set_title(f'SR of asset classes as function of signal for "{signal_type}"')

        if saveout:
            fname = __code_dir__ + f'Analysis/results/port-matrix-using-signal-{signal_type}.pdf'
            print('Saving figure to',fname)
            plt.savefig(fname,bbox_inches='tight')

        return res
        
    ########## replicate stuff in paper to make sure data/code are working ##########
    
    def get_forecast_and_oosR2(used,wt,baseline=None):
        '''
        Use blending weight to construct the blended signal and then calls the SSE-based
        OOS R2.

        baseline -- If this is not None, the use this instead of the rolling mean in the
                    OOS R2 calculation.
        '''
        
        blend = wt * used.pred + (1-wt) * used['mean']

        if baseline is None:
            use_mean = used['mean']
        else:
            use_mean = baseline

        return 100 - 100*((blend-used.true)**2).sum()/((use_mean-used.true)**2).sum()

    def calc_selection_oosR2(self,type,var,show_output=True):
        '''
        Replicate the in-sample model selection charts from Figure 8, the one with
        rows corresponding to each dependent variable, i.e., FutRet, DSpot, DOilVol,
        etc., and columns corresponding to the text, nontext, and both models.
        '''
        
        thed = self.read_lasso_file(type,var)
        
        oosR2s = {}
        eps = 0.01
        wts = np.arange(0,1+eps,eps)

        for lookback in [3,3.5,4,4.5,5]:

            oosR2s_look = {}
            for terc in [1,2,3]:
            
                oosR2 = pd.Series(0.,index=wts)
                for wt in wts:
                    subd = thed[(thed.index <= self.select_end) &
                                (thed[f'lookback_tercile_lag_{lookback:.1f}yr']==terc)]
                    oosR2[wt] = PortEngine.get_forecast_and_oosR2(subd,wt)
                                                                  #baseline=100 if var == 'FutRet' else 0)
                oosR2s_look[terc] = oosR2

            ## collect the data
            oosR2s[lookback] = pd.DataFrame(oosR2s_look)
                
        if show_output:
            ax = oosR2s[4].plot(title=f'Forecasting model selection for {var} and {type}')
            ax.set_ylim(max(-5,oosR2s[4].min().min()),oosR2s[4].max().max()*1.02)
            ax.grid(color='lightgrey',alpha=0.5)

        return thed, oosR2s

    def selection_and_evaluation_oosR2(self,type='text'):
        '''
        Reproduce the text part of Table 6 in the original FAJ submission and extend
        this to allow for the zero (instead of mean) benchmark.
        '''

        rows = []
        opt_idxs = []
        for var in ['FutRet','DSpot','DOilVol','xomRet','bpRet','rdsaRet','DInv','DProd']:

            ## get the three tercile OOSR2 curves (as a function of weight) for each of the
            ## lookbacks
            thed, oosR2s = self.calc_selection_oosR2(type,var,show_output=False)

            row, opt_idx = {'Var':var}, {'Var':var}
            
            for lookback, val in oosR2s.items():

                ## get the maximum R2 points
                max_wt = val.idxmax() ## this gives back a weight for earch tercile (column)
                max_terc = max_wt.idxmax()
                max_wt = max_wt[max_terc]

                row[lookback] = val.loc[max_wt,max_terc]
                opt_idx[lookback] = (max_terc,max_wt)
                
            rows.append(row)
            opt_idxs.append(opt_idx)
            
        rows, opt_idxs = pd.DataFrame(rows), pd.DataFrame(opt_idxs)
        rows.set_index('Var',inplace=True)
        opt_idxs.set_index('Var',inplace=True)
        
        ## get optimal lookback for each row
        rows['L'] = rows.idxmax(axis=1)
        for idx, row in rows.iterrows():
            rows.loc[idx,'Phi'] = opt_idxs.loc[idx,row['L']][0]
            rows.loc[idx,'Wt'] = opt_idxs.loc[idx,row['L']][1]

        ## Get the evaluation set out-of-sample R2s, as in Table 6. For each dependent
        ## variable, e.g., FutRet, DSpot, etc., and a given lookback window, tercile level,
        ## and weight, find the OOS R2s (which is completely out-of-sample). Only keep the
        ## observation if the lagged tercile was in the right bucket.
        print(f'\nGetting evaluation window OOS R2s >= {self.eval_start}:')
        for var in rows.index:

            thed = self.read_lasso_file(type,var)
        
            ## The "true" return starts on the date indicated by index and the pred and mean
            ## use information prior to and including this date.
            used = thed[(thed.index >= self.eval_start) &
                        (thed[f'lookback_tercile_lag_{rows.loc[var,"L"]:.1f}yr'] == rows.loc[var,'Phi'])]
            oosR2 = PortEngine.get_forecast_and_oosR2(used,rows.loc[var,'Wt'])
            rows.loc[var,'OOS R2'] = oosR2

        print(rows.round(2))
        return rows
        
    def check_forward_and_mean_returns(self,var):
        '''
        Sanity check on true returns.
        '''

        thet = self.read_lasso_file('text',var)
        then = self.read_lasso_file('nontext',var)
        theb = self.read_lasso_file('both',var)

        ## make sure all three files are the same
        assert (thet.true == then.true).all()
        assert (thet.true == theb.true).all()

        ## calculate fwd returns to double check
        if var not in ['FutRet','xomRet', 'rdsaRet', 'bpRet']: return
        level = 100 if var == 'FutRet' else 0 ## adjust levels for FutRet only

        ## get the fwd returns to double check
        for tt in thet.index:
            rets_fwd1d = self.rets_fwd1d[var].loc[(tt+pd.Timedelta(4,'w')):(tt+pd.Timedelta(12,'w'))]
            thet.loc[tt,'check fwd'] = (1+rets_fwd1d/100).prod()*100 - (100 if var != 'FutRet' else 0)

        ## get the historical five-year returns to double check
        ## adjust for: 8 week returns; FutRet level; mean level does not have dividends
        ## but the series I am using from BBG contain dividends
        dvdyld = 6
        thet['check mean'] = self.mktd[var].rolling(248*5).mean()*40 \
            + (100 if var == 'FutRet' else 0) \
            - (dvdyld*40/260 if var != 'FutRet' else 0)
        
        ##
        ##  plotting
        ##
        fig, axs = plt.subplots(2,1,figsize=(7,5))
        thet[['true','check fwd']].plot(title=f'{var}: Forward 8-week ahead returns',ax=axs[0])
        thet[['mean','check mean']].plot(title=f'{var}: Backward 5-year mean of 8-week returns' + \
                                         f' assumed {dvdyld}% divd',ax=axs[1])
        plt.tight_layout()
