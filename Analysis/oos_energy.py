import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

__data_loc__ = os.getenv('HOME')+'/code/Energy/OutOfSample'

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
        
    def calc_actual_and_verify(self):
        '''
        This calculates the actual outcomes series as the prediction minus the difference,
        and then verifies that the implied actual series from the different model versions
        (const, base, text, full) are the same within tolerance.
        '''

        print('\tcalc\'ing and verifying actual series')
        
        ## get the implied actual series (not that for a given depser, this should be the same
        ## across all types (i.e., const, base, text, and full forecasting methods
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
                    ii += 1


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
        
    def __repr__(self):

        return f'pdata: price-based timing data {self.pdata.shape}\n' + \
            f'pred8: predictions from 8-wk 1 non-txt & 1 txt model {self.pred8.shape}\n'

    def blended_oos(self,numvar='1_1',mtype='full'):
        '''
        mtype -- which model to compare to const: base, text, full
        '''

        sers = []
        for depser in self.pred8.index.get_level_values('depser').unique():
        
            const = self.pred8.xs((numvar,depser,'const'),level=['numvar','depser','type'])['diff']
            model = self.pred8.xs((numvar,depser,mtype),level=['numvar','depser','type'])['diff']
            ##zero = self.pred8.xs((numvar,depser,'const'),level=['numvar','depser','type'])['actual']
        
            ## align after droppoing na's
            good_idx = const.notna() & model.notna()
            const = const[good_idx]
            model = model[good_idx]
            
            ## get the R2's or RMSE ratios
            evalstat = []
            wts = np.arange(0,1,0.01)
            for wt in wts:
                blend = wt * model + (1-wt) * const
                ##evalstat.append(np.sqrt((blend**2).mean()/(const**2).mean()))
                evalstat.append(100-100*(blend**2).mean()/(const**2).mean())

            sers.append(pd.Series(evalstat,index=wts,name=depser))

        ##
        ##  Plot only up to a max_wt range
        ##
        sers = pd.concat(sers,axis=1)
        fig = plt.figure()
        max_wt = 0.5
        axs = sers[sers.index <= max_wt].plot(subplots=True,layout=(2,4),figsize=(11,5),
                                              xlabel='$w$',ylabel='$R^2_{OOS}$ in percent')
        for ii, ax in enumerate(axs.flatten()):
            ##ax.axhline(1,linestyle='--')
            ax.axhline(0,linestyle='--')
            ax.set_title(f'Global opt = {sers.iloc[:,ii].max():.2f}%')

        model_name = {'base':'non-text','text':'text','full':'full'}
        plt.suptitle(f'Comparing {model_name[mtype]} ({numvar}) and const model blends',y=0.97)
        plt.tight_layout()
        
