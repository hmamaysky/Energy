import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
__code_dir__ = f'c:/users/{os.getenv("USERNAME")}/code/energy/'

class PortEngine:

    def __init__(self):
        ## from Kaiwen Hou email on 5/16/2025:
        ## 1) Is 0-1 the text model and 1-0 the nontext model?
        ## Yes. 
        ## In https://github.com/hmamaysky/Energy/blob/master/OutOfSample/train_demo.ipynb:
        ## df.to_csv(f"res_Forward_Lasso/{num_nontext_var}-{num_text_var}_{d_var}.csv")
        self.fname_map = {'text':'0-1','nontext':'1-0','both':'1-1'}

        ## other params to match the original FAJ submission OOS analysis
        self.select_end = pd.Timestamp('2010-01-01') - pd.DateOffset(weeks=8)
        self.eval_start = '2010-01-01'

        ## Shell PLC trades on Euronext Amsterdam; BP is the ADR return; FutRet is return on
        ## USO ETF; XOM is from NYSE (the data are from Bloomberg). These are daily returns
        ## including dividends. And gb03 is the 3-month T-bill rate. All data are in percent.
        fname = __code_dir__+'Analysis/support/bbg_energy_with_divd_rets 2025-11-19.csv'
        print('Reading',fname)
        self.mktd = pd.read_csv(fname,parse_dates=['Dates'],index_col=0)
        self.dpy = 365.25/self.mktd.index.to_series().diff().dt.days.mean()
        assert 258 < self.dpy < 262 ## assert returns are daily

        ## convert DSpot from a level to a percent return series
        self.mktd.DSpot = self.mktd.DSpot.pct_change(1)*100
        
        ## line up with fwd returns & double check timing <-- these come from the MKTD data
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
        thed = pd.read_csv(fname,parse_dates=True,index_col=0)

        ## get the historical five-year returns to double check
        ## adjust for: 8 week (40 biz day) returns; FutRet level; mean level does not
        ## have dividends but the series I am using from BBG contain dividends
        self.dvdyld = 6
        eight_wks = 40
        thed['mean_check'] = self.mktd[var].rolling(248*5).mean()*eight_wks \
            - (self.dvdyld*eight_wks/260 if var not in ['FutRet','DSpot'] else 0)

        ## add momentum measures
        one_month = 21
        def cagr(xx,lag=None): return (1+xx[:lag]/100).prod()*100-100
        thed['cagr_l1']    = self.mktd[var].rolling(one_month* 1).apply(cagr,raw=True)
        thed['cagr_l3']    = self.mktd[var].rolling(one_month* 3).apply(cagr,raw=True)
        thed['cagr_l6']    = self.mktd[var].rolling(one_month* 6).apply(cagr,raw=True)
        thed['cagr_l12']   = self.mktd[var].rolling(one_month*12).apply(cagr,raw=True)
        thed['cagr_l12m1'] = self.mktd[var].rolling(one_month*12).apply(lambda xx: cagr(xx,-one_month),raw=True)

        ## make futures adjustment
        if var == 'FutRet':
            for el in thed.columns:
                if el[:5] in ['mean_','cagr_']:
                    thed[el] += 100
        
        return thed

    def __repr__(self):

        repstr = ''
        
        for el in ['fname_map','select_end','eval_start','dpy','dvdyld']:
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
        assert signal_type in ['pred','true','mean','cagr_l1','cagr_l3','cagr_l6','cagr_l12','cagr_l12m1',
                               'text_and_cagr']
        assert retser in self.mktd.columns
        print(f'\nPortfolio simulation for {retser} using {type}/{signal_type} model for {signal_var}')

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
        if signal_type == 'text_and_cagr':
            portd['signal'] = 0.9*thed['pred'] + 0.1*thed['cagr_l12m1'] - level
        else:
            portd['signal'] = thed[signal_type] - level

        ## set the weights from the selected signals
        for ii,tt in enumerate(portd.index):

            ## if date exists in the signal data, then set the weights from it
            if not np.isnan(portd.signal[tt]):
                if portd.signal[tt] > +thresh:
                    portd.loc[tt,'weights'] = +1
                elif portd.signal[tt] < -thresh:
                    portd.loc[tt,'weights'] = -1

            ## wights persist from the last signal
            if np.isnan(portd.loc[tt,'weights']):
                if ii > 0:
                    portd.loc[tt,'weights'] = portd.weights.iloc[ii-1]
                else:
                    portd.loc[tt,'weights'] = 1 ## if not signal, start at weight of 1

        ## calculate each day's portfolio return
        assert not portd.weights.isna().any()
        port_rets = portd.weights * portd[retser]
        port_rets.name = 'port_rets'

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
            portd[['signal','weights']].clip(upper=2,lower=-2).ffill().plot(ax=axs[1],xlabel='')
            portd['signal'].clip(upper=2,lower=-2).plot(ax=axs[1],xlabel='',marker='.',color='C0')
            axs[1].axhline(0,color='lightgrey',linestyle='--')
            num_years = (portd.weights.index[-1] - portd.weights.index[0]).days / 365.25
            num_trades = portd.weights.diff()[portd.weights.diff() != 0].shape[0]
            axs[1].set_title(f'Number of trades per year = {num_trades/num_years:.1f}')
            plt.tight_layout()

            if saveout:
                fname = __code_dir__ + f'Analysis/results/port-test-{type}-{signal_var}-{retser}.pdf'
                print('Saving figure to',fname)
                plt.savefig(fname,bbox_inches='tight')
            
        return thed, pd.concat([portd.weights,port_rets],axis=1), sr(portd[retser]), sr(port_rets)

    def port_test_matrix(self,type,signal_type,saveout):
        '''
        type -- which signal to use: text, nontext, both (this is not used for signal_type=='mean')
        signal_type -- the signal (mean, pred, true -- forward looking) to use in trade simulation
        saveout -- save the output?
        '''

        assert type in ['text','nontext','both']
        
        forecast_vars = ['FutRet','bpRet','rdsaRet','xomRet','DSpot']
        dep_vars = ['FutRet','bpRet','rdsaRet','xomRet']

        ## get results
        res = pd.DataFrame(np.nan,index=forecast_vars,columns=dep_vars)
        resl = pd.Series(np.nan,index=dep_vars)
        for forecast in forecast_vars:
            for dep in dep_vars:
                _, _, srl, sr = self.port_test(type,forecast,dep,thresh=0,
                                               signal_type=signal_type,show_plots=False)
                resl[dep] = srl
                res.loc[forecast,dep] = sr

        ## display results
        res.loc['Long'] = resl
        print(res)

        ## plot the return matrix
        plt.figure()
        ax = sns.heatmap(res,annot=True,fmt='.3f',vmin=-0.2,vmax=0.8,cmap='RdYlGn',annot_kws={'size':12})
        ax.set_ylabel('Signal')
        ax.set_xlabel('Return series')
        signal_str = f'"{signal_type}" signal' if signal_type[:4] in ['mean','cagr'] \
            else f'"{type}" model signal'
        ax.set_title(f'SR of asset classes as function of {signal_str}')

        if saveout:
            fname = __code_dir__ + f'Analysis/results/port-matrix-using-signal-{type}-{signal_type}.pdf'
            print('Saving figure to',fname)
            plt.savefig(fname,bbox_inches='tight')

        return res

    def port_test_signals(self,type,signal_var,retser,saveout):
        '''
        type -- one of 'text','nontext','both'
        signal_var -- the variable from which to get the signal
        retser -- the variable to trade
        saveout -- save output?
        '''
        
        allres = {}
        for sigtype in ['pred','cagr_l1','cagr_l3','cagr_l6','cagr_l12','cagr_l12m1','text_and_cagr']:
            ld,res,_,sr = self.port_test('text',signal_var,retser,thresh=0,signal_type=sigtype,show_plots=False)
            allres[sigtype+f' SR: {sr:.3f}']= res.port_rets
        allres = pd.DataFrame(allres)

        ## plotting
        ax = (1+allres/100).cumprod().dropna().plot(title=f'Signals from {signal_var} to {retser}',
                                                    xlabel='',figsize=(9,6))
        ax.legend(ncol=2)
        if saveout:
            addstr = f'{type}-{signal_var}-for-{retser}'
            fname = __code_dir__ + f'Analysis/results/port-test-signals-{addstr}.pdf'
            print('\nSaving figure to',fname)
            plt.savefig(fname,bbox_inches='tight')

        ## for correlation plot, rename cols to dtop SRs
        allres.columns = [el.split(' SR:')[0] for el in allres.columns]
            
        ## correlations
        plt.figure()
        ax = sns.heatmap(allres.corr(),annot=True,fmt='.3f')
        ax.set_title(f'Correlation of return signals from {signal_var} for {retser}')
        ax.text(-0.20, -0.30, f'Data from {allres.dropna().index[0].date()} '+ \
                f'to {allres.dropna().index[-1].date()}', transform=ax.transAxes)

        if saveout:
            fname = __code_dir__ + f'Analysis/results/port-test-signals-{addstr}-corrs.pdf'
            print('Saving figure to',fname,'\n')
            plt.savefig(fname,bbox_inches='tight')

        
    ########## replicate stuff in paper to make sure data/code are working ##########

    def get_benchmark(oos_benchmark,var):

        assert oos_benchmark in ['rolling','zero']

        if oos_benchmark == 'rolling':
            return None
        else:
            return 100 if var == 'FutRet' else 0
    
    def get_forecast_and_oosR2(used,wt,benchmark=None):
        '''
        Use blending weight to construct the blended signal and then calls the SSE-based
        OOS R2.

        benchmark -- If this is not None, the use this instead of the rolling mean in the
                     OOS R2 calculation.
        '''

        ## determine the benchmark (i.e., either the rolling mean, or the passed-in level)
        if benchmark is None:
            benchmark = used['mean']

        ## get the OOS R2
        blend = wt * used.pred + (1-wt) * used['mean']
        return 100 - 100*((blend-used.true)**2).sum()/((benchmark-used.true)**2).sum()

    def calc_selection_oosR2(self,type,var,oos_benchmark='rolling',show_output=True):
        '''
        Replicate the in-sample model selection charts from Figure 8, the one with
        rows corresponding to each dependent variable, i.e., FutRet, DSpot, DOilVol,
        etc., and columns corresponding to the text, nontext, and both models.

        oos_benchmark -- one of 'rolling' (rolling mean in lasso estimation window) or 'zero'
                         (return/change forecast is assumed to be zero)
        '''

        ## determine the appropriate benchmark
        benchmark = PortEngine.get_benchmark(oos_benchmark,var)
        
        ## run the analysys
        thed = self.read_lasso_file(type,var)
        
        oosR2s = {}
        wts = np.linspace(0,1,101)

        for lookback in [3,3.5,4,4.5,5]:

            oosR2s_look = {}
            for terc in [1,2,3]:
            
                oosR2 = pd.Series(0.,index=wts)
                for wt in wts:
                    subd = thed[(thed.index <= self.select_end) &
                                (thed[f'lookback_tercile_lag_{lookback:.1f}yr']==terc)]
                    oosR2[wt] = PortEngine.get_forecast_and_oosR2(subd,wt,benchmark)
                oosR2s_look[terc] = oosR2

            ## collect the data
            oosR2s[lookback] = pd.DataFrame(oosR2s_look)

        ## show output for the 4-year lookback for terciles
        if show_output:

            ## label the series for each tercile
            locR2 = oosR2s[4].copy()
            wt_idx = locR2.idxmax()
            locR2.columns = [f'$\phi$={el} R2={locR2.loc[wt_idx[el],el]:.1f} w={wt_idx[el]:.2f}'
                             for el in locR2.columns]
            
            ## plotting            
            ax = locR2.plot(title=f'Forecasting model selection for {var} and {type}')
            ax.set_ylim(max(-5,oosR2s[4].min().min()),oosR2s[4].max().max()*1.02)
            ax.grid(color='lightgrey',alpha=0.5)

        return thed, oosR2s

    def selection_and_evaluation_oosR2(self,type='text',oos_benchmark='rolling'):
        '''
        Reproduce the text part of Table 6 in the original FAJ submission and extend
        this to allow for the zero (instead of mean) benchmark.

        oos_benchmark -- one of 'rolling' (rolling mean in lasso estimation window) or 'zero'
                         (return/change forecast is assumed to be zero)
        '''

        rows = []
        opt_idxs = []
        for var in ['FutRet','DSpot','DOilVol','xomRet','bpRet','rdsaRet','DInv','DProd']:

            ## get the three tercile OOSR2 curves (as a function of weight) for each of the
            ## lookbacks
            thed, oosR2s = self.calc_selection_oosR2(type,var,oos_benchmark,show_output=False)

            row, opt_idx = {'Var':var}, {'Var':var}
            
            for lookback, val in oosR2s.items():
            
                ## get the maximum R2 points
                max_wt = val.idxmax() ## this gives back the max weight for earch tercile (column)
                max_oosR2 = val.max() ## this is the maximum OOS R2
                max_terc = max_oosR2.idxmax()
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
            benchmark = PortEngine.get_benchmark(oos_benchmark,var)

            ## The "true" return starts on the date indicated by index and the pred and mean
            ## use information prior to and including this date.
            used = thed[(thed.index >= self.eval_start) &
                        (thed[f'lookback_tercile_lag_{rows.loc[var,"L"]:.1f}yr'] == rows.loc[var,'Phi'])]
            oosR2 = PortEngine.get_forecast_and_oosR2(used,rows.loc[var,'Wt'],benchmark)
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

        ## manually check the cagrs
        print('Checking cagr calculations:',end=' ')
        one_month = 21
        for nm in thet.columns:
            if nm in ['cagr_l1','cagr_l3','cagr_l6','cagr_l12','cagr_l12m1']:
                print(nm,end=' ')
                idx_test = self.mktd.index.tolist().index(thet.index[-1])
                if nm[-4:] != '12m1':
                    lag = int(nm.replace('cagr_l',''))
                    skip = 0
                else:
                    lag = 12
                    skip = 21

                ## the +1 is to make inclusive
                start = idx_test-int(one_month*lag)+1
                end = idx_test - skip + 1
                test = (1+self.mktd[var].iloc[start:end]/100).prod()*100
                if var != 'FutRet':
                    test -= 100

                ## check values match
                assert abs(thet[nm].iloc[-1]-test) < 1e-15
        print()
                
        ## calculate fwd returns to double check
        if var not in ['FutRet','DSpot','xomRet','rdsaRet','bpRet']:
            print(f'Variable [{var}] not recognized. Returning.')
            return
        level = 100 if var == 'FutRet' else 0 ## adjust levels for FutRet only

        ## get the fwd returns to double check
        for tt in thet.index:
            rets_fwd1d = self.rets_fwd1d[var].loc[(tt+pd.Timedelta(4,'w')):(tt+pd.Timedelta(12,'w'))]
            thet.loc[tt,'fwd_check'] = (1+rets_fwd1d/100).prod()*100 - (100 if var != 'FutRet' else 0)

        ##
        ##  plotting
        ##
        fig, axs = plt.subplots(3,1,figsize=(9,7.5))
        thet[['true','fwd_check']].plot(title=f'{var}: Forward 8-week ahead returns',ax=axs[0])
        thet[['mean','mean_check']].plot(title=f'{var}: Backward 5-year mean of 8-week returns' + \
                                         (f' assume {self.dvdyld}% divd' if var not in ['FutRet','DSpot'] else ''),
                                         ax=axs[1])
        thet[['cagr_l1','cagr_l3','cagr_l6','cagr_l12','cagr_l12m1']].plot(ax=axs[2])
        plt.tight_layout()

        return thet
