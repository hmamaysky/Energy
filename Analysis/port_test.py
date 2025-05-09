import numpy as np
import pandas as pd
__work_dir__ = '~/code/energy/outofsample/res_Forward_lasso/'

class PortEngine:

    def __init__(self):

        self.fname_map = {'text':'0-1','nontext':'1-0','both':'1-1'}
        
    def read_file(self,type,var):

        assert type in self.fname_map.keys()
        
        ## read the data
        fname = __work_dir__ + self.fname_map[type] + '_' + var + '.csv'
        print('Reading',fname)
        return pd.read_csv(fname,parse_dates=True,index_col=0)

    def __repr__(self):

        repstr = f'fname_map: {self.fname_map}\n'
        
        for el in ['thed']:
            if hasattr(self,el):
                repstr = f'thed: data of size {getattr(self,el).shape}\n'

        return repstr
        
    def test_oosR2(self,type,var):

        thed = self.read_file(type,var)
        
        oosR2s = {}
        wts = np.arange(0,1.05,0.05)

        for terc in [1,2,3]:

            oosR2 = pd.Series(0.,index=wts)
            
            for wt in wts:
        
                subd = thed[(thed.index <= '2009-11-30') &
                            (thed['lookback_tercile_lag_4.0yr']==terc)]
                blend = wt * subd.pred + (1-wt) * subd['mean']
                oosR2[wt] = 100 - 100*((blend-subd.true)**2).sum()/((subd['mean']-subd.true)**2).sum()

            oosR2s[terc] = oosR2

        oosR2s = pd.DataFrame(oosR2s)
        ax = oosR2s.plot()
        ax.set_ylim(max(-5,oosR2s.min().min()),oosR2s.max().max()*1.02)
        ax.grid(color='lightgrey',alpha=0.5)
