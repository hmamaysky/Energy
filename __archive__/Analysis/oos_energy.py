import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm

__data_loc__ = os.getenv('HOME')+'/code/Energy/OutOfSample'
__out_loc__ = os.getenv('HOME')+'/code/Energy/Analysis/results'

class OOSAnalysis:

    def __init__(self):

        print('OOS R2 analysis...')

        fname = __data_loc__ + '/blended/data/transformed_data_prices_v19.dta'
        print('\treading',fname)
        self.pdata = pd.read_stata(fname)

        ## read in the prediction and error files
        def read_preds_errs_file(dtype,numvar):

            fname = f'{__data_loc__}/blended/results/Lasso_{dtype}_10fold_8wk_{numvar}.xlsx'
            print('\treading',fname)
            rawd = pd.read_excel(fname)

            ## convert all cells in file to ordered series
            rawd.set_index('Unnamed: 0',inplace=True)
            sers = []
            len_ser = None
            ## the rows are called either "const_pred" or "const_diff", etc.
            for rr in ['const','base','text','full']:
                for cc in rawd.columns: ## these are the series being forecasted
                    ## parse out the data from the text string
                    preds = [float(el) for el in
                             rawd.loc[rr+'_'+dtype,cc].replace('[','').replace(']','').split(',')]
                    preds = np.array(preds)

                    if cc == 'FutRet' and dtype == 'pred':
                        preds -= 100 ## for some reason these are based at 100
                    
                    ## make sure length of all forecasts is the same
                    if len_ser is None:
                        len_ser = len(preds)
                    else:
                        assert len_ser == len(preds)

                    ## store the data
                    sers.append(pd.DataFrame({'numvar':numvar,'type':rr,'depser':cc,
                                              'index':range(len(preds)),dtype:preds}))

            ## concatenate all the series together
            print(f'\t\tlength of all series = {len_ser}')
            return pd.concat(sers,axis=0)

        diff8 = pd.concat([read_preds_errs_file('diff','1_1'),
                           read_preds_errs_file('diff','2_2')],axis=0)
        pred8 = pd.concat([read_preds_errs_file('pred','1_1'),
                           read_preds_errs_file('pred','2_2')],axis=0)
        self.pred8 = diff8.merge(pred8,on=['numvar','type','depser','index'])

        ## now get the implied actual series
        self.calc_actual_and_verify()
        
    def __repr__(self):

        return f'pdata: price-based timing data {self.pdata.shape}\n' + \
            f'pred8: predictions from 8-wk 1-1 (base/txt/both) and 2-2 (base/txt/both) models {self.pred8.shape}\n'

    def calc_actual_and_verify(self):
        '''
        This calculates the actual outcomes series as the prediction minus the difference,
        and then verifies that the implied actual series from the different model versions
        (const, base, text, full) are the same within tolerance.
        '''

        print('\tcalc\'ing and verifying actual series')
        
        ## Get the implied actual series. Note that for a given depser, this should be the same
        ## across all types (i.e., const, base, text, and full forecasting methods.
        self.pred8['actual'] = self.pred8['pred'] - self.pred8['diff']
        self.pred8.set_index(['index','numvar','type','depser'],inplace=True)

        ## check that this is the correct calculation for actual series for all
        ## (depser,numvar,type) tuples, e.g., ('FutRet','2_2','const')
        for depser in self.pred8.index.get_level_values('depser').unique():
            ser = None
            ii = 0
            ## for a given depser, all actual series should be the same
            for mtype in self.pred8.index.get_level_values('type').unique():
                for numvar in self.pred8.index.get_level_values('numvar').unique():
                    if ii == 0:
                        act1 = self.pred8.xs((depser,mtype,numvar),
                                             level=['depser','type','numvar'])['actual'].to_list()
                    else:
                        act2 = self.pred8.xs((depser,mtype,numvar),
                                             level=['depser','type','numvar'])['actual'].to_list()
                        diff = np.array([abs(el1 - el2) for el1,el2 in zip(act1,act2)
                                         if (not np.isnan(el1) and not np.isnan(el2))])
                        if diff.max() > 1e-13:
                            print('There is a problem reconstructing actual series for',
                                  (depser,mtype,numvar),f'. Max diff is {diff.max()}.')
                            raise Exception('Reconstructed actual series do not reconcile.')
                    ii += 1

            print(f'\t\t {depser:7}: implied "actual" series good for ' + \
                  f'{{const,text,base,full}} x {{1_1,2_2}}')
                    
    def gen_old_table(self,numvar='1_1'):
        '''
        Generate the old Table VI from the paper, which shows the RMSE ratios.
        '''

        def calc_vr(ser,null):
            return np.sqrt((ser**2).mean()/(null**2).mean())
        
        rows = []
        for el in self.pred8.index.get_level_values('depser').unique():

            used = self.pred8.xs((numvar,el),level=['numvar','depser'])
            const = used.xs('const',level='type')['diff']
            base = used.xs('base',level='type')['diff']
            text = used.xs('text',level='type')['diff']
            full = used.xs('full',level='type')['diff']

            ## we can use the actual outcome series compared against the zero model,
            ## i.e., the forecast is always zero, as an alternative to the const (i.e.,
            ## rolling mean) model
            zero = used.xs('const',level='type')['actual']
            
            rows.append({'depser':el,
                         #'zero_const':calc_vr(const,zero),
                         #'zero_base':calc_vr(base,zero),
                         #'zero_text':calc_vr(text,zero),
                         #'zero_full':calc_vr(full,zero),
                         'const_text':calc_vr(text,const),
                         'const_full':calc_vr(full,const),
                         'non-text_text':calc_vr(text,base),
                         'non-text_full':calc_vr(full,base)})

        table = pd.DataFrame(rows)
        print(f'Table for {numvar}:')
        print(table.round(3),end='\n\n')
        return table            
        
    def check_const_for_actual(self):
        '''
        Check how well constant forecasts actual outcomes across the dependent
        variables.
        '''

        numvar = '1_1' ## should be same for '2_2'
        for depser in self.pred8.index.get_level_values('depser').unique():

            ## calc the regression of actual on constant model
            actual = self.pred8.xs((numvar,depser,'const'),level=['numvar','depser','type'])['actual']
            c_pred = self.pred8.xs((numvar,depser,'const'),level=['numvar','depser','type'])['pred']
            mod = sm.OLS(actual,sm.add_constant(c_pred),missing='drop')
            res = mod.fit()
            print(f'depser = {depser:8s} (numvar={numvar})  R2={res.rsquared:6.4f}' + \
                  f'  bb={res.params["pred"]:6.3f}  pval={res.pvalues["pred"]:6.3f}')
            

    def blended_oos(self,numvar='1_1',mtest='full',mcontrol='const',saveout=False):
        '''
        mtest -- which model to compare to const: base, text, full
        mcontrol -- which is the control model, i.e., const for main results, but can also be
                    base if we want to use the non-text model as the control model for OOS R2
        '''

        oos_R2 = []
        mean_err = []
        for depser in self.pred8.index.get_level_values('depser').unique():
        
            ## these are the errors of the control model and the rolling lasso model
            model = self.pred8.xs((numvar,depser,mtest),level=['numvar','depser','type'])['diff']
            control = self.pred8.xs((numvar,depser,mcontrol),level=['numvar','depser','type'])['diff']

            ## this is the prediction coming from the control model, e.g., rolling mean for const
            control_pred = self.pred8.xs((numvar,depser,mcontrol),level=['numvar','depser','type'])['pred']

            ## this is the actual outcome for the given depser
            actual = self.pred8.xs((numvar,depser,mtest),level=['numvar','depser','type'])['actual']

            ## this is if we want to set a fixed mean forecast
            ##use_mean = -1

            ## align after droppoing na's
            good_idx = control.notna() & model.notna() & control_pred.notna() & actual.notna()
            model = model[good_idx]
            control = control[good_idx]
            control_pred = control_pred[good_idx]
            actual = actual[good_idx]

            ## sanity check (this is redundant given 'calc_actual_and_verify'
            assert np.abs(control - (control_pred - actual)).max() < 1e-13
            
            ## get the R2's or RMSE ratios
            oos_R2_local = []
            mean_err_local = []
            wts = np.arange(0,1,0.01)
            for wt in wts:
                blend = wt * model + (1-wt) * (control_pred-actual)

                oos_R2_local.append(100-100*(blend**2).mean()/((control_pred-actual)**2).mean())

                mean_err_local.append(blend.mean()/actual.std())

            ## collect results for the given depser
            oos_R2.append(pd.Series(oos_R2_local,index=wts,name=depser))
            mean_err.append(pd.Series(mean_err_local,index=wts,name=depser))

        ## collect all the data
        oos_R2 = pd.concat(oos_R2,axis=1)
        mean_err = pd.concat(mean_err,axis=1)

        ##
        ##  Plot only up to a max_wt range
        ##
        def plot_outcome(data,label,ylabel,hline=None):

            fig = plt.figure()
            if mcontrol == 'const':
                max_wt = 0.50
            else:
                max_wt = 1
            axs = data[data.index <= max_wt].plot(subplots=True,layout=(2,4),figsize=(11,4),
                                                  xlabel='$w$',ylabel=ylabel)
            if hline is not None:
                for ii, ax in enumerate(axs.flatten()):
                    ax.axhline(hline,linestyle='--')
                    ax.set_title(f'Global opt = {data.iloc[:,ii].max():.2f}%')

            model_name = {'base':'non-text','text':'text','full':'full','const':'rolling mean'}
            title = f'{label} {model_name[mtest]} ({numvar}) vs {model_name[mcontrol]} model blends'
            plt.suptitle(title,y=0.97)
            plt.tight_layout()

            if saveout:
                fname = (__out_loc__ + '/' + title + '.pdf').\
                    replace('(','').replace(')','').replace(' ','-').replace('_','-')
                print('Saving to',fname)
                plt.savefig(fname)

        plot_outcome(oos_R2,'OOS R2',ylabel='$R^2_{OOS}$ in percent',hline=0)
        plot_outcome(mean_err,'Mean normalized error (pred-actual)',ylabel='Normalized mean error')
