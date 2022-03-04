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
# %% read in and process OOS results
oos = en.OOSResults()
# %% calculate various stats
oosall = oos.calc('All',saveout=True)
#oostxt = oos.calc('Text',saveout=True)
# %% check the successful model distributions
oos.variable_distribution(1)
aa = en.correlated_binom(45,0,0.5,common=0.5,nsims=250)
# %% count the run types
def foo(res,commons=np.arange(0,1.1,0.1)):
    
    ex = {}
    for common in commons:
        
        subres = res[res.index.str.contains(f'sim{common:.1f}$')]
        ex[f'sim{common:.1f}'] = \
            subres[(subres >= '(0.95)') | (subres <= '(0.05)')].notna().sum().sum()

    ex = pd.Series(ex)
    print(ex)
    
    return ex
    
# %% simulations to verify closed form probs of length-k runs
df = oos.compare_sim_data(num_sims=50000)
# %% probs for runs -- example
prun = oos.prob_of_run(0.282,3,7,verbose=True)
print(prun)
# %% some data plots
en.plot_depvar_corrs()