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
oos.calc()
# %% probs for runs -- example
from scipy.stats import binom
qq = 0.339
kk = 4
nrun = 16
maxruns = 45*22
prun = oos.prob_of_run(0.282,kk,7)
print(qq,'Prob run of length {} = {}'.format(kk,prun))
nruns = range(max(0,nrun-100),min(nrun+100,maxruns))
plt.plot(nruns,binom.pmf(nruns,maxruns,prun))
print('Prob > {} runs = {}'.format(nrun,\
      1-sum(binom.pmf(range(nrun+1),maxruns,prun))))

############################## Run the simulation for oos subperiod analysis and create a comparison table ##############################

import energy as en
oos = en.OOSResults()
df = oos.compare_sim_data() # m=100000 is set as default
df1 = oos.compare_sim_data(m=1000000)
## save the results as csv files.
df.to_csv('sim_result.csv')
df1.to_csv('sim_result1.csv')


