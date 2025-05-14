import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
__code_dir__ = '~/code/energy/'

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

        ## line up with fwd returns & double check timing
        self.rets_fwd1d = self.mktd.shift(-1)
        assert self.rets_fwd1d.xomRet['2025-05-08'] == self.mktd.xomRet['2025-05-09']
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

    def port_test(self,type,var,retser,lag,tercile,thresh=2,show_plots=True):
        '''
        type -- which signal to use: text, nontext, both
        var -- which variable's lasso model to use, e.g. FutRet, DSpot, etc.
        retser -- which market to trade, one of FutRet, xomRet, rdsaRet, and bpRet
        lag -- which lag tercile to use
        tercile -- value tercile at lag should have to trade
        thresh -- threshold of forecast for taking a position
        show_plots -- do the visualization or not
        '''

        ## sanity checks
        assert retser in self.mktd.columns
        
        ## get the lasso signals file
        thed = self.read_lasso_file(type,var)
        level = 100 if var == 'FutRet' else 0 ## levels across series
        usesd = thed.true.std()
        assert (level-0.2*usesd) <= thed.true.mean() <= (level+0.2*usesd) ## make sure have the correct level
        
        ## get the excess returns and then calculate the trading weights
        ##idx = (self.rets_fwd1d.index >= self.eval_start) & (self.rets_fwd1d.index <= thed.index[-1])
        idx = (self.rets_fwd1d.index >= thed.index[0]) & (self.rets_fwd1d.index <= thed.index[-1])
        portd = pd.DataFrame(self.rets_fwd1d[retser][idx] - self.rets_fwd1d.gb03[idx],columns=[retser])
        portd['weights'] = np.nan

        ## get the signal
        ## used = thed[thed[f'lookback_tercile_lag_{lag:.1f}yr'] == tercile]
        portd['signal'] = thed['pred'] - level

        ## print status
        weight_AR = 0.
        print(f'Portfolio simulation for {retser} using {type} model for {var} AR={weight_AR}')

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
                    portd.loc[tt,'weights'] = 0

            #if ii == 0:
            #    portd.loc[tt,'weights'] = change_to_weight
            #else:
            #    portd.loc[tt,'weights'] = weight_AR * portd.weights.iloc[ii-1] + change_to_weight

        ## calculate each day's portfolio return
        assert not portd.weights.isna().any()
        port_rets = portd.weights * portd[retser]

        sr = port_rets.mean()/port_rets.std()*np.sqrt(self.dpy)
        print(f'Port rets SR = {sr}')

        if show_plots:
        
            ## get the underlying return series
            jnt_rets = pd.DataFrame({retser:portd[retser],'port':port_rets})
            cum_rets = (1+jnt_rets/100).cumprod()
            cum_rets.columns = [f'{el}: {cum_rets[el].iloc[-1]:.2f} ' + \
                                f'SR: {jnt_rets[el].mean()/jnt_rets[el].std()*np.sqrt(self.dpy):.3f}'
                                for el in jnt_rets.columns]

            fig, axs = plt.subplots(2,1,figsize=(8,6))
            cum_rets.plot(title=f'Cumulative excess returns for {retser} using {var}\n' + \
                          f'signal {type} with lag={lag} and tercile={tercile}',
                          ax=axs[0],xlabel='')
            portd[['weights','signal']].ffill().clip(upper=2,lower=-2).plot(ax=axs[1],xlabel='')
            axs[1].axhline(0,color='lightgrey',linestyle='--')
        
        return thed, portd.weights, sr
        
    ########## replicate stuff in paper to make sure data/code are working ##########
    
    def get_forecast_and_oosR2(used,wt):
        '''
        Use blending weight to construct the blended signal and then calls the SSE-based
        OOS R2.
        '''
        
        blend = wt * used.pred + (1-wt) * used['mean']
        return 100 - 100*((blend-used.true)**2).sum()/((used['mean']-used.true)**2).sum()

    def select_oosR2(self,type,var):
        '''
        Replicate the in-sample model selection charts from Figure 8, the one with
        rows corresponding to each dependent variable, i.e., FutRet, DSpot, DOilVol,
        etc., and columns corresponding to the text, nontext, and both models.
        '''
        
        thed = self.read_lasso_file(type,var)
        
        oosR2s = {}
        wts = np.arange(0,1.05,0.05)

        for terc in [1,2,3]:

            oosR2 = pd.Series(0.,index=wts)
            
            for wt in wts:
        
                subd = thed[(thed.index <= self.select_end) &
                            (thed['lookback_tercile_lag_4.0yr']==terc)]
                oosR2[wt] = PortEngine.get_forecast_and_oosR2(subd,wt)

            oosR2s[terc] = oosR2

        oosR2s = pd.DataFrame(oosR2s)
        ax = oosR2s.plot(title=f'Forecasting model selection for {var} and {type}')
        ax.set_ylim(max(-5,oosR2s.min().min()),oosR2s.max().max()*1.02)
        ax.grid(color='lightgrey',alpha=0.5)

    def eval_oosR2(self,type,lag,tercile,weight,var):
        '''
        Get the evaluation set out-of-sample R2s, as in Table 6. For each dependent
        variable, e.g., FutRet, DSpot, etc., and a given lookback window, tercile level,
        and weight, find the OOS R2s (which is completely out-of-sample). Only keep the
        observation if the lagged tercile was in the right bucket.
        '''

        thed = self.read_lasso_file(type,var)

        ## XXXHM: what is the timing of this date? Does the "true" return start on this date?
        used = thed[(thed.index >= self.eval_start) &
                    (thed[f'lookback_tercile_lag_{lag:.1f}yr'] == tercile)]
        oosR2 = PortEngine.get_forecast_and_oosR2(used,weight)
        print(f'Evaluation OOS R2 for {var} {type} = {oosR2:.4f}\n')
        
        return used

    def check_true_rets(self,var):
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

        ## get the fwd returns to double check
        for tt in thet.index:
            rets_fwd1d = self.rets_fwd1d[var].loc[(tt+pd.Timedelta(4,'w')):(tt+pd.Timedelta(12,'w'))]
            thet.loc[tt,'Check'] = (1+rets_fwd1d/100).prod()*100

        ## adjust the level (only FutRet shows gross returns)
        if var != 'FutRet':
            thet['Check'] -= 100
            
        thet[['true','Check']].plot(title=f'{var}')

        
