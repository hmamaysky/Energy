from Energy.Analysis import energy as en
# %% read in text data
txtd = en.read_info()
# %% read in the processed data
serd = en.read_processed()
en.summary(serd)  ## summary of series
# %% search for outlier event
en.find_outliers(serd,txtd)
# %% display outlier events
en.show_outliers(serd,txtd)
# %% some text series plots
en.plot_text_stats(serd)
# %% check SR for futret series (can run not on grid)
import os
locdir = os.getenv('HOME')+'/code/energy/data/'
vard = pd.read_stata(locdir+'transformed_data_prices_v19.dta')
fr = vard.FutRet_t8_Fri/100-1
print('Mean = {:.4f}  SD = {:.4f}  SR = {:.4f}'.\
      format(fr.mean(),fr.std(),np.sqrt(52/8)*fr.mean()/fr.std()))
# %% *** can run locally, not on grid ***
## read in and process OOS results about number of runsin subperiods
oos = en.OOSResults()
# %% sample run (for debugging run_runs_sims.py)
res = oos.res_for_depvar('bpRet','Text')
# %% summarize outputs of [run_runs_sims.py]
oosall = oos.summary('All',saveout=True)
oostxt = oos.summary('Text',saveout=True)
# %% check empirical successful variable distributions
oos.variable_distribution(runlen=1)
# %% simulated successful variable distributions
numA = len(oos.model_summary()['all_vars'])
numB = 0
pr = 0.5 ## probability of a successful trial
en.correlated_binom(numA,numB,pr,common=0.,nsims=2500)
en.correlated_binom(numA,numB,pr,common=0.9,nsims=2500)
# %% simulations to verify closed form probs of length-k runs
df = oos.compare_sim_data(num_sims=50000)
# %% probs for runs -- example
prun = oos.prob_of_run(0.282,3,7,verbose=True)
print(prun)
# %% some data plots
en.plot_depvar_corrs()