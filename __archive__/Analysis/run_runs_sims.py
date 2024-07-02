#!/c/apps/anaconda3/python

from Energy.Analysis import energy as en
import re, os, datetime, multiprocessing

__out_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/results/'

def run_one_depvar(oos,dv,varset,run_date):

    ## do the actual calculation
    res = oos.res_for_depvar(dv,varset)

    ## save to csv
    fname = __out_dir__ + f'sims-for-runs-tests-{varset}-{dv}-{run_date}.csv'
    print('Saving to',fname)
    res.T.to_csv(fname)
    
def run_varset(oos,varset,run_date):

    print(f'Running for {varset}')

    ## kick off runs
    pool = multiprocessing.Pool(processes=len(oos.data.depvar.unique()))
    for dv in oos.data.depvar.unique():
        pool.apply_async(run_one_depvar,(oos,dv,varset,run_date))
    pool.close()
    pool.join()

    return(1)
    
if __name__ == '__main__':

    print('Running simulations to calculate runs p-values with correlations...')

    oos = en.OOSResults()
    
    ## run for either 'All' or 'Text' forecasting variables
    run_date = str(datetime.date.today())
    run_varset(oos,'All',run_date)
    run_varset(oos,'Text',run_date)
    
