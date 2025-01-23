import pandas as pd
import os, re, numpy as np
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.stats import binom, norm
from scipy.linalg import cholesky
from collections import Counter
from utils import pyhelp as pyh
import statsmodels.api as sm

__text_dir__ = '/shared/share_mamaysky-glasserman/energy_drivers/2020-11-16'
__out_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/results'
__sup_dir__ = os.getenv('HOME')+'/code/Energy/Analysis/support'
__dpi_res__ = 500


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
        
    def variable_distribution(self,runlen):
        '''
        Get the empirical distribution of forecasting variables appearing in models
        with runs. This may help with calibration of the correlation coefficient in
        the correlated binomial simulations.

        runlen -- Length of run to look for variable distribution
        '''

        runcol = f'Run{runlen}'
        assert runcol in self.data.columns

        ## summary stats about models
        summ = self.model_summary()
        numA = len(summ['all_vars'])
        numT = len(summ['txt_vars'])

        for model_type in ['all','text']:

            if model_type == 'all':
                num_models = numA * (numA-1) / 2
                var_set = summ['all_vars']
            else:
                num_models = numT * (numT-1) / 2 + numT * (numA - numT)
                var_set = summ['txt_vars']
                
            ## go through all the dependent variables
            hists = []
            for dv in sorted(set(self.data.depvar)):

                ## get all data for given dependent variable with rows w/out missing observations,
                ## and make sure to restrict attention to the desired set of forecasting variables
                dd = self.data[(self.data.depvar==dv)&(self.data[runcol])&
                               (self.data.var1.apply(lambda el: el in var_set) |
                                self.data.var2.apply(lambda el: el in var_set))]

                vs = dd['var1'].to_list() + dd['var2'].to_list()
                vs = pd.Series(Counter(vs)).sort_values()
                vs.name = f'{dv}  std={vs.std():.2f}'

                hists.append(pd.Series(vs).sort_values())

            ## combine hists
            hists = pd.concat(hists,axis=1)
            hists.hist(figsize=(10,7),rwidth=0.8)
            title = f'Distribution of number of successful {model_type}-variable models'

            ## some info about model
            title2 = f'Run length={runlen}  Num vars={numA}  Txt vars={numT}  Num models={num_models:.0f}'

            ## make title and save
            plt.suptitle(title + '\n' + title2,y=0.99,fontsize=15)

            fname = __out_dir__+f'/Dist-{model_type}-'+title2.replace(' ','-')+'-'+str(date.today())+'.png'
            print('Saving to',fname)
            plt.savefig(fname,dpi=__dpi_res__)
        

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

    def res_for_depvar(self,dv,varset):
        '''
        Run the analysis for a single dependent variable.
        '''

        ## number of simulations and common variation for correlated outcomes case
        nsims = 100
        nsims = 2500
        commons = np.arange(0,1.1,0.2)
        commons_str = ', '.join([str(el.round(2)) for el in commons])
        print(f'For {dv} using c (common) values of {commons_str}.')
        
        ## get vars
        thed = self.data[(self.data.depvar==dv) & (self.data.drop_row != True)].copy(deep=True)
        summ = self.model_summary()
        txt_vars = summ['txt_vars']
        all_vars = summ['all_vars']
        nontxt_vars = summ['nontxt_vars']
        run_cols = summ['run_cols']
        
        ## status update
        print('Working on:',dv,end='\n\n')

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
               '# All/Txt':'{}/{}'.format(len(all_vars),len(txt_vars)),
               'nsims':nsims,
               'commons':'-'.join([f'{el:.3f}' for el in commons])}

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

            ## run simulations to get spread (i.e., measure of dispersion of frequency counts
            ## of number of times variables appear in two-variable forecasting models) and
            ## simulated p-values
            spread_means = []  
            pval_sims = []
            for common in commons:
                sims, spreads = correlated_binom(varsA,varsB,prun,common=common,nsims=nsims,plot=False)
                pval_sim = len(sims[sims > nrun])/len(sims)

                spread_means.append(spreads.mean())
                pval_sims.append(pval_sim)

            ## compare spread
            if rr != 'Run0':
                subd = thed[thed[rr]>0]
                succvars = subd['var1'].tolist() + subd['var2'].tolist()
            else: ## this is the case of 0-length runs

                ## get all pairs of variables for models
                varpairs = []
                if varset == 'All':
                    for ii in range(len(all_vars)):
                        for jj in range(ii+1,len(all_vars)):
                            varpairs += [(all_vars[ii],all_vars[jj])]
                else: ## working with Text variable
                    for ii in range(len(txt_vars)):
                        for jj in range(ii+1,len(txt_vars)):
                            varpairs += [(txt_vars[ii],txt_vars[jj])]
                    for ii in range(len(txt_vars)):
                        for jj in range(len(nontxt_vars)):
                            varpairs += [(txt_vars[ii],nontxt_vars[jj])]

                ## drop those pairs that have a run
                for kk in range(thed.shape[0]):

                    flag1 = (thed.iloc[kk].var1,thed.iloc[kk].var2) in varpairs
                    flag2 = (thed.iloc[kk].var2,thed.iloc[kk].var1) in varpairs

                    if flag1:
                        varpairs.remove((thed.iloc[kk].var1,thed.iloc[kk].var2))

                    if flag2:
                        varpairs.remove((thed.iloc[kk].var2,thed.iloc[kk].var1))

                ## count the occurrences of the remaining pairs that had runs of
                ## length 0
                succvars = [el[0] for el in varpairs] + [el[1] for el in varpairs]

            ## now the code proceeds the same for either Run0 or non-Run0 runs
            succvars = Counter(succvars)

            for el in all_vars:
                if el not in succvars.keys():
                    succvars[el] = 0

            ## sanity check
            assert sum(succvars.values())/2 == nrun

            ## proximity to simulations
            spread = np.std(list(succvars.values()))
            idx = np.argmin(abs(np.array(spread_means)-spread))

            print(f'{rr}, best common {commons[idx]}, pval {pval_sims[idx]}')

            res[rr+'-pval-sim'] = f'({pval_sims[idx]:.2f})'
            res[rr+'-common'] = f'[{commons[idx]:.1f}]'

        res = pd.DataFrame(res,index=[0]).set_index('Dep Var')
        return res

    def model_summary(self):
        '''
        Some summary statistics about the OOS model tests.
        '''

        ## the columns containing the run counts
        run_cols = [el for el in self.data.columns if re.search('Run[0-9]',el)]
        
        ## get the vars which aren't missing in any subperiod
        alld = self.data[self.data.drop_row != True]
        txt_vars = set(alld[alld.var1_text==1].var1).union(alld[alld.var2_text==1].var2)
        all_vars = set(alld.var1).union(alld.var2)
        nontxt_vars = all_vars - txt_vars
        
        ## convert to  lists
        txt_vars, all_vars, nontxt_vars = list(txt_vars), list(all_vars), list(nontxt_vars)
        
        ## By focusing on drop_row == False, there should be no missing values anywhere.
        all_vars_str = '|'.join(all_vars)
        print(f'Checking that there are no missing observations among {len(all_vars)} variables: ',end='')
        num_missing = self.data[(self.data.drop_row == True) &
                                self.data.var1.str.contains(all_vars_str) &
                                self.data.var2.str.contains(all_vars_str)].shape[0]
        assert num_missing == 0
        print('All good!\n')

        ## return vals
        return {'all_vars':all_vars,
                'txt_vars':txt_vars,
                'nontxt_vars':nontxt_vars,
                'run_cols':run_cols}
    
    def summary(self,varset,saveout=False):
        '''
        Calculate the statistics for the simulations to determine runs p-values.

        varset -- Either 'All' to look at results for all the forecasting variables or 'Text'
        to look at the results for models that include at least one text variables.
        saveout -- Save/plot the table form of the results?
        '''

        assert varset in ['All','Text']
        
        ## get all results
        fnames = pyh.most_recent_file(__out_dir__,f'sims-for-runs-tests-{varset}',multiple=True)
        allres = []
        for fname in fnames:
            allres.append(pd.read_csv(fname,index_col=0))
        allres = pd.concat(allres,axis=1)
            
        ## same some results that will be discarded when displaying
        nsims = allres.loc['nsims'].unique()
        assert len(nsims)==1
        nsims = nsims[0]
        
        commons = allres.loc['commons'].unique()
        assert len(commons)==1
        commons = [float(el) for el in commons[0].split('-')]

        ## reorder rows and add all rows whose name starts with "Run[0-9]"
        rows = ['Runs','Max Runs','q','# All/Txt']

        ## find strings that end with "Run[0-9]+"
        runs = allres.index[allres.index.str.contains('Run[0-9]+$')]
        for run in runs:
            add_rows = allres.index[allres.index.str.contains(run)]
            rows.extend(add_rows)
        
        ## reorder how dependent variables show up
        allres = allres[['FutRet','DSpot','DOilVol','xomRet','bpRet','rdsaRet','DInv','DProd']].loc[rows]

        ## make columns names into first row
        toprow = pd.DataFrame(dict(zip(allres.columns,allres.columns)),index=['Dep Var'])
        allres = pd.concat([toprow,allres],axis=0)
        
        ## show/save output?
        try:
            if saveout:
                self.plot_calc(allres,varset,nsims,commons,detailed=False)
        except:
            print('Something went wrong with plot_calc.')
            
        return allres

    def plot_calc(self,allres,varset,nsims,commons,detailed=False):
        '''
        Plot the results of the calc() method.
        '''
        
        assert varset in ['All','Text']

        ## select which rows
        if not detailed:
            allres = allres[allres.index.str.contains('common') == False]
            ##allres = allres[allres.index.str.contains('pval$') == False]

        ## convert to plotting
        pltres = allres.reset_index()

        ## prettify row labels
        pltres['index'] = [re.sub('Run[0-9]-','',el) for el in pltres['index']]
            
        fig = plt.figure(figsize=(7,0.25*pltres.shape[0]),dpi=__dpi_res__)
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
        fname = __out_dir__ + \
            f'/runs-tests-for-{varset}-sims-{nsims}-comm-{len(commons)}-{date.today()}.png'
        print('Saving to',fname)
        plt.savefig(fname,bbox_inches='tight',dpi=__dpi_res__)
    
    def compare_sim_data(self, str_len=7, num_sims=100000):
        """
        Compares the probability of seeing at least one run of length kk with the fraction of
        rows that have at least one run of length kk calculated by simulation.
        """
        calc_res = self.summary(varset='All')
        qlist = [float(el) for el in calc_res.loc['q']]
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
        sim_df.loc['q', :] = qlist
        
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
    num_models = int(varsA*(varsA-1)/2 + varsA*varsB)
    print(f'Using: cutoff = {cutoff:.3}  common = {common:.2}  # models = {num_models}')
    
    ## run sims
    res = []
    spread = []
    for aa in range(nsims):

        if (aa+1) % 500 == 0:
            print(aa,end=' ')
        
        ## draw bunch of correlated normals
        rv = norm.rvs(size=(num_models,1))
        fA = norm.rvs(size=(varsA,1))
        
        ## create factor structure of A-type shocks
        base = 0
        idx1 = []  ## store the index of 1st variable go into the outcome
        idx2 = []  ## store the index of 2nd variable go into the outcome
        for ii in range(varsA-1):
            locf = np.sqrt(0.5) * (fA[ii] + fA[(ii+1):])
            ref = range(base,base + varsA - ii - 1)
            rv[ref] = np.sqrt(common) * locf + np.sqrt(1-common) * rv[ref]
            base = base + varsA - ii - 1

            ## which indexes do these represent
            idx1.extend([f'Var{ii}']*len(ref))
            idx2.extend([f'Var{el}' for el in range(ii+1,varsA)])

        ## add to these the B-type shocks if there are any
        if varsB != 0:
            fB = norm.rvs(size=(varsB,1))
            for ii in range(varsA):
                locf = np.sqrt(0.5) * (fA[ii] + fB) ## one A-type and all B-types
                ref = range(base,base + varsB)
                rv[ref] = np.sqrt(common) * locf + np.sqrt(1-common) * rv[ref]
                base = base + varsB

                ## which indexes do these represent
                idx1.extend([f'Var{ii}']*varsB)
                idx2.extend([f'Var{el}' for el in range(varsA,varsA+varsB)])

        ## sanity check
        assert base == num_models ## the count variables (base) and ex-ante calc'd (num_models)
        assert len(idx1) == num_models
        assert len(idx2) == num_models
                
        ##
        ##  Number of successes: check number less than cutoff
        ##
        n_success = len(rv[rv<cutoff])
        res.append(n_success)

        ##
        ##  Distribution of successful forecating pairs
        ##
        hists = [idx1[ii] for ii, xx in enumerate(rv < cutoff) if xx] \
            + [idx2[ii] for ii, xx in enumerate(rv < cutoff) if xx]
        hists = Counter(hists)

        ## not sure this is necessary? maybe?
        for ii in range(varsA+varsB):
            vname = f'Var{ii}'
            if vname not in hists.keys():
                hists[vname] = 0

        ## save the number of successes and a statistic about the spread
        spread.append(np.std(list(hists.values())))

    print()
        
    ## plotting
    if plot:

        run_date = date.today()
        
        ## plot the histogram of successful models
        plt.figure()
        plt.hist(res,rwidth=0.9,density=True,bins=50)
        xxr = range(np.min(res),np.max(res)+1)
        zz = binom.pmf(xxr,num_models,prob_success)
        plt.plot(xxr,zz,color='red')
        plt.suptitle(f'Number successful trials with common={common:.2f} nsims={nsims}',y=0.94)
        
        fname = __out_dir__ + f'/number-successful-trials-common-{common:.2f}-{date.today()}.png'
        print(fname)
        plt.savefig(fname,dpi=__dpi_res__)
        
        ## plot the histogram of spreads between occurrences of most and least
        ## frequently successful models
        plt.figure()
        plt.hist(spread,rwidth=0.9,color='navy',
                 label=f'common = {common:.2f}\nnsims = {nsims}\nmean = {np.mean(spread):.2f}')
        plt.suptitle('Spread in occurrence between high and low frequency models',y=0.94)
        plt.legend(handlelength=0,handletextpad=0)
        
        fname = __out_dir__ + f'/spread-high-low-frequency-{common:.2f}-{date.today()}.png'
        print(fname)
        plt.savefig(fname,dpi=__dpi_res__)

        
    return np.array(res), np.array(spread)

    
############################## Read in text data ##############################

def assert_on_grid():

    if not os.getenv('USER'):
        raise Exception('This has to be run on the computing grid, not on PC.')
        

def read_info():

    assert_on_grid()
    
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

    assert_on_grid()
    
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
    

############################## the regression analysis ##############################

class EnergyBetas:

    def __init__(self,sdate='2000-01-01'):

        self.read_bbg(sdate)
        self.read_french(sdate)
        
    def read_french(self,sdate):

        f5 = pyh.FrenchReader.read_five_factor()
        mom = pyh.FrenchReader.read_momentum()

        ## drop early dates
        f5 = f5[f5.index >= sdate]
        mom = mom[mom.index >= sdate]
        momf = mom[['PRIOR 8','PRIOR 9','Hi PRIOR']].mean(axis=1) - \
            mom[['Lo PRIOR','PRIOR 2','PRIOR 3']].mean(axis=1)

        ## agg to monthly
        f5 = f5.resample('ME').apply(lambda xx: (1+xx/100).prod()*100 - 100)
        momf = momf.resample('ME').apply(lambda xx: (1+xx/100).prod()*100 - 100)
        momf.name = 'UMD'
        
        ## combine data
        self.factors = pd.concat([f5,momf],axis=1)

    def read_bbg(self,sdate):

        self.sdate = sdate
        self.rets_map = {'XLE':'Energy',
                         'XLB':'Basics/materials',
                         'XLI':'Industrials',
                         'XLC':'Communications',
                         'XLP':'Cons staples',
                         'XLF':'Financials',
                         'XLV':'Healthcare',
                         'XLK':'Technology',
                         'XLRE':'REITs',
                         'XLU':'Utilities',
                         'XLY':'Cons discretionary',
                         #'SPY':'S&P 500',
                         #'IEF':'7-10-yr Treasury',
                         'KBE':'Banks',
                         'IWM':'Russell 2000',
                         'VXUS':'Int\'l ex-US',
                         'VGK':'Europe',
                         'EEM':'Emerging Markets',
                         #'MDY':'S&P MidCap 400',
                         'VTV':'Value',
                         'MTUM':'Momentum'}
        
        ## get the returns data
        print('BBG: Getting returns data...')
        tiks = [(el+' equity','DAY_TO_DAY_TOT_RETURN_GROSS_DVDS',el) for el in self.rets_map.keys()]
        self.rets_mo = pyh.bdh_data(tiks,start_date=self.sdate.replace('-',''),
                                    periodicity='MONTHLY')

        ## get the levels data, including generic 3-month T-bill
        print('BBG: Getting levels data...')
        levels = [('CL1 comdty','PX_LAST','Oil price'),('GB1M index','PX_LAST','RF')]
        sdate_levels = (pd.Timestamp(sdate) - pd.Timedelta(days=32)).date()
        self.levels = pyh.bdh_data(levels,start_date=str(sdate_levels).replace('-',''),
                                   periodicity='MONTHLY')

    def calc(self):

        print('Calculting some relationships...')
        
        ## add to returns data
        self.rets_mo['Oil'] = self.levels['Oil price'].pct_change(1)*100
        self.rets_mo['RF'] = self.levels.RF.shift(1)/12

        print('Full sample annualized vols:')
        print(self.rets_mo.std()*np.sqrt(12))

    def __repr__(self):

        ret_str = f'sdate: {self.sdate}\n'
        for el in ['rets_map','rets_mo','levels','factors']:
            if hasattr(self,el):
                ret_str += f'{el} {type(getattr(self,el))}: has {len(getattr(self,el))} entries/rows\n'

        return ret_str
        
    def plot(self):
        '''
        Show cumulative returns of FF6 factors.
        '''
        
        (1+self.factors/100).cumprod().plot(subplots=True,figsize=(8,6),layout=(4,2))

    def regs(self):
        '''
        Run regressions to show oil betas.
        '''

        print('Running OLS factor regressions...')

        ## get the factors (the FF ones and oil)
        factors = pd.concat([self.factors[[el for el in self.factors.columns if el != 'RF']],
                             self.rets_mo.Oil],axis=1)

        factors = factors.apply(lambda xx: xx * factors['Mkt-RF'].std() / xx.std(),axis=0)
        
        ## add constant
        factors['const'] = 1
        print(factors.std())
        
        ## run all regressions
        allres = {}
        for tik in [el for el in self.rets_mo.columns if el not in ['Oil','RF']]:

            ## get excess returns and drop NAs
            exret = self.rets_mo[tik] - self.factors.RF
            alld = pd.concat([factors,exret],axis=1).dropna()
            
            ## get the Xs and Ys
            XX = alld.iloc[:,:-1]
            yy = alld.iloc[:,-1]

            ## run OLS with Newey-West with lags = int(TT^0.25);
            ##   >> for the functional form, see Hoechle, 2007, "Robust Standard Errors," Stata
            n_lags = int(4*(len(yy)/100)**(2/9))
            print(f'\t{tik}: N-W lags={n_lags} start={XX.index[0].date()}...')

            ## model w/ oil
            model = sm.OLS(yy,XX)
            res = model.fit(cov_type='HAC',cov_kwds={'maxlags':n_lags})

            ## model w/out oil
            model_wo = sm.OLS(yy,XX[[el for el in XX.columns if el != 'Oil']])
            res_wo = model_wo.fit()
           
            ## save results
            get_star = lambda pv: '***' if pv <= 0.01 else '**' if pv <= 0.05 else '*' if pv <= 0.1 else ''
            store_params = [f'{coeff:.3f}{get_star(pval)}' for coeff,pval in zip(res.params,res.pvalues)]
            allres[tik] = pd.Series(store_params,index=res.params.index)
            allres[tik]['R2'] = f'{res.rsquared:.3f}'
            allres[tik]['Del R2'] = (res.rsquared - res_wo.rsquared).round(3)

        ## combine results & drop some rows
        allres = pd.DataFrame.from_dict(allres,orient='index')
        allres = allres[[el for el in allres.columns if el not in ['const','Del R2']]]
        allres = allres.sort_values('Oil',ascending=False,
                                    key=lambda xx: xx.str.replace('*','').astype(float))

        ## rename indexes & reduce to values
        allres.index = [self.rets_map[el] for el in allres.index]
        allvals = allres.apply(lambda xx: xx.str.replace('*','').astype(float),axis=1)
        print(allvals)

        ## plot
        plt.figure(figsize=(9.5,7))
        ax = sns.heatmap(allvals,annot=allres,fmt='s')
        ax.xaxis.tick_top()
        ax.set_title('Traded securities: factor and oil betas',y=1.06,fontsize=14)
        ax.text(0,-0.05,transform=ax.transAxes,fontsize=12,
                s=f'Monthly data {self.rets_mo.index[0].date()} to {self.rets_mo.index[-1].date()}' + \
                ' | Factor vols normalized to equal Mkt-Rf')
