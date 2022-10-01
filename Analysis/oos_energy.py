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
                
    def __repr__(self):

        return f'pdata: price-based timing data {self.pdata.shape}\n' + \
            f'pred8: predictions from 8-wk 1 non-txt & 1 txt model {self.pred8.shape}\n'

    def dopred(self,numvar='1_1',depser='FutRet'):

        used = self.pred8[['numvar','depser','type','diff']]
        const = used[(self.pred8.numvar == numvar) &
                     (self.pred8.depser == depser) &
                     (self.pred8.type == 'const')]['diff']
        full = used[(self.pred8.numvar == numvar) &
                    (self.pred8.depser == depser) &
                    (self.pred8.type == 'full')]['diff']

        ############################## align after droppoing na's

        ## get the RMSE ratios
        ratios = []
        wts = np.arange(0,0.25,0.01)
        for wt in wts:

            blend = pd.Series(wt * full.to_numpy() + (1-wt) * const.to_numpy())
            ratios.append(np.sqrt((blend**2).mean()/(const**2).mean()))

        print(ratios)
        plt.plot(wts,ratios)
