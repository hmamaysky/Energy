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
oosatt = oos.calc('All',saveout=True)
oostxt = oos.calc('Text',saveout=True)
# %% simulations to verify closed form probs of length-k runs
df = oos.compare_sim_data(num_sims=500000)
# %% probs for runs -- example
prun = oos.prob_of_run(0.282,3,7,verbose=True)
print(prun)
# %% some data plots
en.plot_depvar_corrs()