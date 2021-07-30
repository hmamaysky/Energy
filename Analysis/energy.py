import pandas as pd
import os, re, numpy as np
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

__text_dir__ = '/shared/share_mamaysky-glasserman/energy_drivers/2020-11-16'
__out_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/results'
__sup_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/support'

'''
Get total count of length k runs.  A run must either be at the right edge and have
a zero after or at the left edge and have a zero before.  Or it can be in the middle
but then it has to have zeros around it.
Prob of observing a run of length k for a given pair (1 or 1128) is
p = (q**4 * (1-q)**1) * 2 (the corse) + (qq**4 * (1-q)**2) * 2 (the interior)
If this is the prob of having length k-run in any pair, the prob of seeing m
runs across 1128 pairs is the binomial pdf (n m) * p**m * (1-p)**(n-m).
We can then tabulate the prob of each m for each observed LHS variable.
'''


############################## parse OOS persistence results ##############################

class OOSResults():

    def __init__(self):
        
        res = []
        for fname in os.listdir(__sup_dir__+'/OOS-subperiod-results'):

            if not '.csv' in fname: continue

            locres = pd.read_csv(__sup_dir__+'/OOS-subperiod-results/'+fname)
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
        print('\tDuplicate models? ',end='')
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
        print()


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
        
        
    def calc(self):
        '''
        Calculate the statistics.
        '''

        ## the columns containig the run counts
        run_cols = [el for el in self.data.columns if re.search('Run[0-9]',el)]
        
        ## stats for a single dependent variable
        print('Working on:',end=' ')
        def res_for_depvar(thed):

            print(thed.depvar[0],end=' ')

            ## get the cases where the row is not dropped
            thed = thed[thed.drop_row != True]

            ## get the vars the make it
            txt_vars = set(thed[thed.var1_text==1].var1).union(thed[thed.var2_text==1].var2)
            all_vars = set(thed.var1).union(thed.var2)

            ## under assumption each variable gets selected at least once, other than
            ## those that get dropped
            ################################################## what set of variables are ever dropped
            tot_mods = len(all_vars)*(len(all_vars)-1)/2
            tot_txt_mods = (len(all_vars)-len(txt_vars)) * len(txt_vars) \
                + len(txt_vars)*(len(txt_vars)-1)/2

            ## calculate the statistics for the current variable
            num_txt = (thed.var1_text | thed.var2_text).sum()
            res = {'Dep Var':thed.depvar[0],
                   'All Runs':thed.shape[0],
                   'Txt Runs':num_txt,
                   ## probability of beating is the total number of beats divided by the total
                   ## number of possible beats (total possible models x number of periods)
                   'Pr beat':round(thed[self.per_cols].to_numpy().sum() / (tot_mods * len(self.per_cols)),3),
                   'Num All':len(all_vars),
                   'Num Txt':len(txt_vars)
            }

            ## number of occurrences of diff length runs
            counts = thed[run_cols].sum()
            res.update(counts)

            return pd.DataFrame(res,index=[0])

        ## get all results
        allres = self.data.groupby('depvar').apply(res_for_depvar).reset_index(drop=True)
        
        ## combine all results
        allres = allres.fillna(0)
        print()
        print(allres)

        return allres
        
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
