import pandas as pd
import os, re, numpy as np
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
import seaborn as sns
import math
from collections import Counter
from scipy.stats import binom, norm
from scipy.linalg import cholesky

__text_dir__ = '/shared/share_mamaysky-glasserman/energy_drivers/2020-11-16'
__out_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/results'
__sup_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/support'


############################## parse OOS persistence results ##############################

class OOSResults():

    def __init__(self):
        
        res = []
        print('Reading in result files:')
        for fname in os.listdir(__sup_dir__+'/OOS-subperiod-results'):

            if not '.csv' in fname: continue

            full_fname = __sup_dir__+'/OOS-subperiod-results/'+fname
            print('\t'+full_fname)
            locres = pd.read_csv(full_fname)
            locres['depvar'] = fname.replace('.csv','')
            res.append(locres)

        ## some cleanups
        colmap = {"('2003-04-25', '2005-03-31')":'Per1',
                  "('2005-04-01', '2007-11-30')":'Per2',
                  "('2007-12-01', '2009-06-30')":'Per3',
                  "('2009-07-01', '2012-02-29')":'Per4',
                  "('2012-03-01', '2014-10-31')":'Per5',
                  "('2014-11-01', '2017-06-30')":'Per6',
                  "('2017-07-01', '2020-01-31')":'Per7',
                  'var_1':'var1', 'var_1_text':'var1_text',
                  'var_2':'var2', 'var_2_text':'var2_text'}
        self.data = pd.concat(res,axis=0).rename(columns=colmap)

        ## find all text vars
        self.per_cols = [el for el in self.data.columns if re.search('Per[0-9]',el)]

        ## replace -1 observations with NAs
        self.data[self.data[self.per_cols]==-1] = np.nan

        ## check data
        self.check()

        ## calculate numbers of runs
        self.calc_runs()

    def check(self):
        '''
        Some sanity checks on the data.
        '''
        
        print('Checking for validity.')
        
        ## Check #1: should not have same variable in model twice
        assert (self.data['var1'] != self.data['var2']).all()

        ## Check #2: see if any pairs exists where {a,b} shows up as {b,a}
        print('\tChecking for duplicate models in: ',end='')
        for depvar in self.data.depvar.unique():

            print(depvar,end=' ')
            thed = self.data[self.data.depvar == depvar]
            dups = []
            combos = (thed.var1+thed.var2).to_list()
            for depvar, ii, jj in zip(thed.depvar,thed.var1,thed.var2):
                if jj+ii in combos:
                    idx1 = combos.index(ii+jj)
                    idx2 = combos.index(jj+ii)
                    dups.append({'Y':depvar,'v1':ii,'v2':jj,'idx1':idx1,'idx2':idx2})

            assert len(dups) == 0
        print('... none found.')


    def calc_runs(self):
        '''
        Calculate how many runs of each type there are for given {depvar,var1,var2}
        combination.  A run is a series of consecutive 1's.
        '''

        print('Calculating numbers of runs.')
        
        ## get names of indicator columns
        def num_runs(rr):

            ## If there are any missing values in the period counts, then one of the
            ## forecasting variables is missing in this period, and we can't use this
            ## for the OOS runs comparisons since not all periods are observed.  Then
            ## ignore this row for future analysis of runs, and send back None
            if rr.isna().any(): return {'drop_row':True}

            ## convert indicators to a string and then find all runs of 1, i.e.
            ## 1+ -- which means all instances of 1 occurring 1 or more times
            strform = ''.join(str(int(el)) for el in rr[self.per_cols])
            all_runs = re.findall('1+',strform)

            ## now convert run matched strings to lengths, and count them
            all_runs = [len(el) for el in all_runs]
            return Counter(all_runs)

        ## calcualte runs in every row, rename and reorder columns
        runs_counts = self.data.apply(num_runs,axis=1,result_type='expand')
        runs_counts.columns = [el if not isinstance(el,int) else 'Run'+str(el)
                               for el in runs_counts.columns]
        runs_counts = runs_counts[runs_counts.columns.sort_values()]

        self.data = pd.concat([self.data,runs_counts],axis=1)
        
    def prob_of_run(self,qq,kk,nn,verbose=False):
        '''
        Params
        ======
        qq - probability of a 1 (and 1-qq is probability of zero)
        kk - run length
        nn - length of string

        Summary
        =======
        This calculates the probability of seeing a run of 1's of length kk in a string of
        0's and 1's of length nn, when the probability of a 1 is qq.
        A run must either be at the right edge and have a zero after or at the left edge and have a
        zero before.  Or it can be in the middle but then it has to have zeros around it.
        Prob of a corner run of length kk is: q**kk * (1-qq)
        Prob of an interior run of length kk is: (1-qq) * qq**kk * (1-qq)
        The only issue is that sometimes later kk-length runs may overlap with prior ones and so
        their probabilities must be adjusted to preclude a prior kk-length run including the
        current one to avoid double counting.  This is the role of the cose labeled (*).

        If this is the prob pr of having a length kk-run in any pair, the prob of seeing m
        runs across (num_vars choose 2) pairs is the binomial pdf: (n m) * pr**m * (1-pr)**(n-m).
        We can then tabulate the prob of each m for each observed LHS variable.
        '''

        ## all failures (0 runs)
        if kk == 0:
            return (1-qq)**nn
        
        ## all successes (run of length n)
        if kk == nn:
            return qq**nn
        
        ## shorter runs can be at the corners or the interior
        prs = np.zeros(nn-kk+1)

        ## recursion for probabilities
        for ii in range(prs.shape[0]):

            ## check if at corner or interior
            if ii == 0 or ii == (prs.shape[0]-1):
                prs[ii] = qq**kk * (1-qq)          ## the starting left-corner 1's sequence
            else:
                prs[ii] = (1-qq) * qq**kk * (1-qq) ## the 0 followed by 1's followed by zero

            ## (*) now decide if needs to modify to exclude overlapping probabilities
            if ii - (kk+1) >= 0:
                ## this sums up all probabilities at lags strictly prior to current position
                ## minus 2, but for pos-2 have to divide by (1-q) to take out the zero that's
                ## already present in that string
                prs[ii] = (1 - prs[0:(ii-kk-1)].sum() - prs[ii-kk-1]/(1-qq)) * prs[ii]

        if verbose:
            print(prs)

        return prs.sum()

            
    def calc(self,varset,saveout=False):
        '''
        Calculate the statistics.

        varset -- Either 'All' to look at results for all the forecasting variables or 'Text'
        to look at the results for models that include at least one text variables.
        saveout -- Save/plot the table form of the results?
        '''

        assert varset in ['All','Text']

        nsims = 1000 ## number of simulations for correlated outcomes case
        
        ## the columns containing the run counts
        run_cols = [el for el in self.data.columns if re.search('Run[0-9]',el)]
        
        ## get the vars which aren't missing in any subperiod
        thed = self.data[self.data.drop_row != True]
        txt_vars = set(thed[thed.var1_text==1].var1).union(thed[thed.var2_text==1].var2)
        all_vars = set(thed.var1).union(thed.var2)

        ## stats for a single dependent variable
        print('Working on:',end=' ')
        def res_for_depvar(dv):

            thed = self.data[(self.data.depvar==dv) & (self.data.drop_row != True)].copy(deep=True)
            print(dv,end=' ')

            ## specialize results to those with at least a text var?
            if varset == 'Text':
                thed = thed[(thed.var1_text==1) | (thed.var2_text==1)]
                tot_mods = (len(all_vars)-len(txt_vars))*len(txt_vars) + len(txt_vars)*(len(txt_vars)-1)/2
            else:
                tot_mods = len(all_vars)*(len(all_vars)-1)/2

            res = {'Dep Var':dv,
                   'Runs':thed.shape[0],
                   'Max Runs':int(tot_mods),
                   ## probability of beating is the total number of beats divided by the total
                   ## number of possible beats (total possible models x number of periods)
                   'q':round(thed[self.per_cols].to_numpy().sum() / (tot_mods * len(self.per_cols)),3),
                   '# All/Txt':'{}/{}'.format(len(all_vars),len(txt_vars))
            }

            ## put in the zero runs row; 'Runs' contains all combinations with at least
            ## one period of outperformance relative to constant
            res['Run0'] = res['Max Runs'] - res['Runs']

            ## get number of pairs {var1,var2} that have a run of a certain length
            runs_indic = thed[run_cols].apply(lambda xx: (xx >= 1).sum(),axis=0)
            res.update(runs_indic)

            ## calculate the p-value for each number
            for rr in ['Run0'] + run_cols:
                qq = res['q']
                kk = int(re.findall('[0-9]',rr)[0])
                nrun = res[rr]
                prun = self.prob_of_run(qq,kk,nn=len(self.per_cols))
                ##print(qq,'Prob run of length {} = {}'.format(kk,prun))
                pval = 1 - sum(binom.pmf(range(nrun+1),tot_mods,prun))

                ##print('Prob > {} runs = {}'.format(nrun,pval))
                res[rr+'-sprob'] = f'{prun:.3f}'
                res[rr+'-pval'] = f'({pval:.2f})'

                ##
                ## Simulations: check distribution allowing for correlated draws
                ##
                if varset == 'Text':
                    varsA = len(txt_vars)
                    varsB = len(all_vars)-len(txt_vars)
                else:
                    varsA = len(all_vars)
                    varsB = 0

                for common in np.arange(0,1.1,0.1):
                    sims = correlated_binom(varsA,varsB,prun,common=common,nsims=nsims,plot=False)
                    pval_sim = len(sims[sims > nrun])/len(sims)
                
                    res[rr+f'-pval-sim{common:.1f}'] = f'({pval_sim:.2f})'

            return pd.DataFrame(res,index=[0])

        ## get all results
        allres = []
        for dv in self.data.depvar.unique():
            allres.append(res_for_depvar(dv).reset_index(drop=True))
        allres = pd.concat(allres,axis=0)
            
        ## reorder columns (these will become the rows of the table post-transpose)
        cols = ['Dep Var','Runs','Max Runs','q','# All/Txt']
        for ii in range(0,len(run_cols)+1):
            add_cols = allres.columns[allres.columns.str.contains(f'Run{ii}')]
            #cols.append('Run'+str(ii))
            #cols.append('Run'+str(ii)+'-sprob')
            #cols.append('Run'+str(ii)+'-pval')
            #cols.append('Run'+str(ii)+'-pval2')
            cols.extend(add_cols)

        allres = allres[cols]

        ## transpose table -- result cols becomes rows and cols are LHS vars
        allres = allres.fillna(0).transpose()
        allres.columns = allres.loc['Dep Var']
        
        ## reorder how dependent variables show up
        allres = allres[['FutRet','DSpot','DOilVol','xomRet','bpRet','rdsaRet','DInv','DProd']]

        print('\n',allres)

        ## show/save output?
        if saveout:

            pltres = allres.reset_index()
            ## prettify row labels
            pltres['index'] = [re.sub('Run[0-9]-','',el) for el in pltres['index']]
            
            fig = plt.figure(figsize=(7,0.25*pltres.shape[0]),dpi=200)
            tbl = plt.table(pltres.round(3).values,
                            loc='center', colWidths=[1.25] + [1]*(pltres.shape[1]-1),
                            bbox=[0,0,1,1])

            ## set edges
            for ii in range(pltres.shape[0]):
                for jj in range(pltres.shape[1]):

                    if jj == 0:
                        fmt = 'LR'
                    elif jj == pltres.shape[1]-1:
                        fmt = 'R'
                    else:
                        fmt = ''
                        
                    if ii == 0:
                        fmt += 'BT'
                    elif re.search('Run[0-9]$',pltres['index'][ii]):
                        fmt += 'T'
                    elif ii == pltres.shape[0]-1:
                        fmt += 'B'
                    tbl[ii,jj].visible_edges = fmt

            plt.axis('tight')
            plt.axis('off')
            plt.title(f'Out-of-sample runs analysis: {varset} models',
                      y=0.99,fontsize=12)

            ## save output
            fname = __out_dir__ + f'/runs-tests-for-{varset}-sims-{nsims}-{date.today()}.pdf'
            print('Saving to',fname)
            plt.savefig(fname,bbox_inches='tight')

        return allres
    
    def compare_sim_data(self, str_len=7, num_sims=100000):
        """
        Compares the probability of seeing at least one run of length kk with the fraction of
        rows that have at least one run of length kk calculated by simulation.
        """
        calc_res = self.calc(varset='All')
        qlist = calc_res.loc['q']
        klist = [1, 2, 3, 4, 5]

        ## Create a table that summarizes the comparison results
        columns = []
        for qq in qlist:

            ## simulate the draws
            arr = np.random.binomial(1, qq, (num_sims, str_len))
            str_arr = arr.astype(str)  ## array of n strings, each of length 1 
            joined_str = [''.join(row) for row in str_arr] ## array of length-n strings

            ## check empirical qq
            emp_qq = arr.sum() / (num_sims*str_len)
            print('{} Sims: Target q = {:.4f}  Empirical q = {:.4f}'.format(num_sims,qq,emp_qq))
            
            ## extract all occurrences of '1's
            run_counts = []
            for row in joined_str:
                ones = [len(el) for el in re.findall('1+', row)]
                run_counts.append(Counter(ones))

            run_counts = pd.DataFrame(run_counts)

            ## get probability of runs of at least a certain length
            sim_probs = run_counts.notna().sum(axis=0) / len(run_counts)
            sim_probs = sim_probs.sort_index() ## return in increasing order

            ## compare the simulated results to the closed-form results
            rows = []
            for kk in klist:
                prun = self.prob_of_run(qq=qq, kk=kk, nn=str_len)
                diff = prun - sim_probs[kk]

                rows.extend([prun, sim_probs[kk], diff])

            columns.append(rows)

        ## collect all the results
        sim_df = pd.DataFrame(columns).transpose()
        
        index_names = []
        for kk in klist:
            index_names.extend(['prun'+str(kk),
                                'frac_sim'+str(kk),
                                'diff'+str(kk)])
            
        sim_df.index = index_names
        sim_df.loc['q', :] = qlist.values
        
        sim_df = sim_df.loc[['q'] + index_names, :]
        sim_df.columns = calc_res.loc['Dep Var', :]
        return sim_df


def correlated_binom(varsA,varsB,prob_success,common=0.3,nsims=2500,plot=True):
    '''
    Compare the PDF of correlated binomials with non-correlated ones to make sure
    PDF is identical.

    varsA -- the number of possible variables that go into 2-variable models
    (each model is assumed to be characterized by a length-7 string of 0s and 1s);
    each type A variable is assumed to enter with every other type A variable

    varsB -- the number of type B variables; type B variables only interact with type A
    variables, not not with other type B variables (e.g., type A are text, and type B are non-text
    would capture all models containing at least one text variable)
    '''

    ## get normal value corresponding to prob_success
    cutoff = norm.ppf(prob_success)

    ## calc the number of models
    num_models = math.comb(varsA,2) + varsA*varsB
    print(f'Using: cutoff = {cutoff:.3}  # models = {num_models}')
    
    ## run sims
    res = []
    for aa in range(nsims):

        if (aa+1) % 500 == 0:
            print(aa,end=' ')
        
        ## draw bunch of correlated normals
        rv = norm.rvs(size=(num_models,1))
        fA = norm.rvs(size=(varsA,1))
        
        ## create factor structure of A-type shocks
        base = 0
        for ii in range(varsA-1):
            locf = np.sqrt(0.5) * (fA[ii] + fA[(ii+1):])
            idx = range(base,base + varsA - ii - 1)
            rv[idx] = np.sqrt(common) * locf + np.sqrt(1-common) * rv[idx]
            base = base + varsA - ii - 1

        ## add to these the B-type shocks if there are any
        if varsB != 0:
            fB = norm.rvs(size=(varsB,1))
            for ii in range(varsA):
                locf = np.sqrt(0.5) * (fA[ii] + fB) ## one A-type and all B-types
                idx = range(base,base + varsB)
                rv[idx] = np.sqrt(common) * locf + np.sqrt(1-common) * rv[idx]
                base = base + varsB
                
        ## check number less than cutoff
        n_success = len(rv[rv<cutoff])

        res.append(n_success)

    print()
    assert base == num_models ## the count variables (base) and ex-ante calc'd (num_models)

    ## plotting
    if plot:
        plt.hist(res,density=True,bins=50)
        xxr = range(np.min(res),np.max(res)+1)
        zz = binom.pmf(xxr,num_models,prob_success)
        plt.plot(xxr,zz,color='red')

    return np.array(res)
        
    
############################## Read in text data ##############################

def read_info():

    ## list all files
    info_dir = __text_dir__+'/DataProcessing/info/'
    fnames = sorted(os.listdir(info_dir))

    ## read all files
    dfs = []

    for nm in fnames:
        fname = info_dir+nm
        print('Reading "',fname,'"',sep='')

        dfs.append(pd.read_csv(fname,index_col='Id'))

    return pd.concat(dfs)


def read_processed(saveout=True):

    ## 12/16/2020: Hongyu says this is the data we're using for the paper
    fname = __text_dir__ + '/data/transformed_data_prices_v14.dta'
    print('Reading in "',fname,'"',sep='')
    procd = pd.read_stata(fname)

    ## mapping from topic numbers to names
    top_short = {'1':'Co','2':'Gom','3':'Env','4':'Epg','5':'Bbl','6':'Rpc','7':'Ep'}
    top_long = {'1':'Company (Co)','2':'Global Oil Markets (Gom)','3':'Environment (Env)',
                '4':'Energy/Power Generation (Epg)','5':'Crude Oil Physical (Bbl)',
                '6':'Refining & Petrochemicals (Rpc)','7':'Exploration & Production (Ep)'}

    ## save?
    if saveout:

################################################## also take , 

        
        ft1 = ['ftopic'+str(el) for el in [1,2,3,4,5,6,7]]
        st1 = ['stopic'+str(el) for el in [1,2,3,4,5,6,7]]
        subd = procd[['date_thurs','Artcount','Entropy',*ft1,*st1,
                      'artcount_4wk','entropy_4wk',
                      *[el+'_4wk' for el in ft1],
                      *[el+'_4wk' for el in st1],
                     'PCAsent','PCAfreq','PCAall']]
        fname = __out_dir__ + '/text-data-' + str(date.today()) + '.csv'
        print('Saving to',fname)
        subd.to_csv(fname)
    
    ## return the data and the topic map
    return {'data':procd,'top_short':top_short,'top_long':top_long}


def summary(serd):
    
    ## look at the correlation matrix of the 7 frequency and topical sentiment series
    f_nms = ['ftopic{}_4wk'.format(el) for el in serd['top_short'].keys()]
    s_nms = ['stopic{}_4wk'.format(el) for el in serd['top_short'].keys()]

    cm = serd['data'][[*f_nms,*s_nms]].corr()
    plt.figure(figsize=(12,12))
    sns.heatmap(cm,annot=True)  ## show numerical values


############################## display series and articles ##############################

## For a given date and topic number, show recent articles that might indicate
## what's going on.
def show_articles(txtd,date_str,top_num,disp_name):

    lookback_days = 28
    topic_thresh = 0.8
    num_articles = 5
    min_length = 100
    min_entropy = 2
    
    date2 = datetime.strptime(date_str,'%Y-%m-%d')
    date1 = date2 - timedelta(days=lookback_days)
                
    ## find all text articles in this topic within the lookback window
    loctxt = txtd[(txtd.TimeStamp_NY >= str(date1)) &
                  (txtd.TimeStamp_NY <= str(date2)) &
                  (txtd.entropy >= min_entropy) &
                  (txtd.total >= min_length) &
                  (txtd['Topic{}'.format(top_num)] > topic_thresh)].sort_values('sentiment')

    ## show output of article headlines
    rows = []

    rows.append({'Sent':None,'Ent':None,'Date':None,
                 'Headline':'{} from {} to {}'.format(disp_name,str(date1)[:10],str(date2)[:10])})

    for jj in range(min(num_articles,loctxt.shape[0])):
        rows.append({'Sent':'{:.3f}'.format(loctxt.iloc[jj]['sentiment']),
                     'Ent':'{:.3f}'.format(loctxt.iloc[jj]['entropy']),
                     'Date':'{}'.format(loctxt.iloc[jj]['TimeStamp_NY'][:10]),
                     'Headline':'{}'.format(loctxt.iloc[jj]['headline'])})

    ## return data frame
    thed = pd.DataFrame(rows)
    print(thed)
    return thed
        

def plot_text_stats(serd):

    labels = {'artcount_4wk':'Article Counts',
              'entropy_4wk':'Entropy',              
              'PCAsent':'First PC of Normalized Topical Sentiment',
              'PCAfreq':'First PC of Normalized Topical Frequency',
              'PCAall':'First PC of All Normalized Text Variables'}
    
    plt.figure(figsize=(10,9))
    nrow = 3
    ncol = 2

    for ii,ser in enumerate(labels.keys()):

        ax = plt.subplot(nrow,ncol,ii+1)
        ax.plot(serd['data'].date_wed,serd['data'][ser],color='blue')
        ax.set_title(labels[ser])

    plt.suptitle('Panel C: Article Counts, Entropy, and Principal Components',fontsize=18)
    plt.tight_layout()

    fname = __out_dir__+'/Figure-2-Panel-C-'+date.today().strftime('%Y-%m-%d')+'.png'
    print('Saving figure to',fname)
    plt.savefig(fname)


def show_outliers(serd,txtd):
    ''' Show a specific set of dates and outliers. '''

    ## these are manually located using find_outliers().
    event_dates = {'1':('2000-09-20','UK fuel protests'),
                   '2':('2002-04-24','Failed Venezuelan coup'),
                   '3':('2015-10-14','Volkswagen emissions scandal'),
                   '4':('2002-02-13','Post-bankruptcy Enron hearings'),
                   '5':('2005-09-21','Hurricane Katrina'),
                   '6':None,
                   '7':('2010-06-02','BP oil spill aftermath')}

    ## dimensions of the picture
    nrow = 2
    ncol = 4

    ## plot the different panels
    for ser in ('ftopic','stopic'):
        
        print(ser)

        ## to store the table
        rows = []

        ## plotting
        plt.figure(ser,figsize=(14,8))
        
        for top_num in serd['top_short'].keys():

            ## get the column name for the {series,topic} pair, e.g., ftopic4_wk
            top_ser = '{}{}_4wk'.format(ser,top_num)
            
            ax = plt.subplot(nrow,ncol,int(top_num))
            ax.plot(serd['data'].date_wed,
                    serd['data'][top_ser],color='blue')
            ax.set_title(serd['top_long'][top_num])

            ## add event points
            if event_dates[top_num] != None:
                event_date = event_dates[top_num][0]
                ax.plot(datetime.strptime(event_date,'%Y-%m-%d'),
                        serd['data'].loc[serd['data'].date_wed == event_date,
                                         top_ser],
                        marker='*',color='red',markersize=16)

                ## get the rows of the headline tables
                rows.append(show_articles(txtd,event_date,top_num,
                                          '{}: {}'.format(serd['top_short'][top_num],
                                                    event_dates[top_num][1])))
                
        ## more plot formatting
        plt.suptitle('Panel A: Topical Frequency' if ser == 'ftopic'
                     else 'Panel B: Topical Sentiment', fontsize=18)
        plt.tight_layout()

        ## save the pics and show the tables
        if ser == 'stopic':

            ## tables
            thed = pd.concat(rows)
            fname = __out_dir__+'/sample-sentences-'+date.today().strftime('%Y-%m-%d')+'.csv'
            print('Saving to',fname)
            thed.to_csv(fname)

            fname = __out_dir__+'/Figure-2-Panel-B-'+date.today().strftime('%Y-%m-%d')+'.png'
            
        else:

            fname = __out_dir__+'/Figure-2-Panel-A-'+date.today().strftime('%Y-%m-%d')+'.png'

        ## pics
        print('Saving figure to',fname)
        plt.savefig(fname)
            
            
############################## find outlier dates ##############################

def find_outliers(serd,txtd):
    ''' Find outlier dates and articles. '''
    
    ## find the sentiment 4wk columns
    ## ntop -- number of top outliers
    ## lookback -- number of days to look back
    def series_outliers(col_str,label,sort_ascending):

        ntop = 2    ## number of top events to examine
        n_diff = 4  ## number of weeks over which to difference for identifying outliers
        
        ## find the columns to search through
        topics = serd['top_short'].keys()

        ## go through all the columns
        for top_num in topics:

            ## get the data corresponding to this series
            top_ser = '{}{}_4wk'.format(col_str,top_num)
            top_name = '{} {}'.format(label,serd['top_short'][top_num])
            used = serd['data'][['date_wed',top_ser]].copy()
            used['dSer'] = used[top_ser].diff(n_diff)
            used = used.sort_values('dSer',ascending=sort_ascending)

            ## keep output string
            print('***',top_name,'***',end='\n\n')
            
            ## go through the more extreme months in each
            for ii in range(ntop):

                ## these are the {news-series,date} pairs for outlier events
                show_articles(txtd,str(used.iloc[ii]['date_wed'])[:10],top_num,top_name)
                
                ## the output string
                print()

    series_outliers('stopic','Sent',sort_ascending=True)

############################## plot some summary stats ##############################

def plot_depvar_corrs():

    dd = pd.read_stata('c:/users/harry/code/Energy/data/transformed_data_physical_v19.dta')
    cc = dd[['DProd_Wed','DInv_Wed','DSpot_Tue','FutRet_Tue',
             'DOilVol_Tue','xomRet_Tue','bpRet_Tue','rdsaRet_Tue']].corr()
    cc.columns = cc.columns.str.replace('_Wed','')
    cc.columns = cc.columns.str.replace('_Tue','')
    cc.index = cc.index.str.replace('_Wed','')
    cc.index = cc.index.str.replace('_Tue','')
    sns.heatmap(cc,annot=True,fmt='.2f')
    ## show the dates of the analysis too
    print('From {} to {}'.format(dd['date_Tue'][0],dd['date_Tue'][1]))
    
