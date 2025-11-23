from Energy.Analysis import port_test as pt
saveout = False
# %% read data
pe = pt.PortEngine()
thed = pe.check_forward_and_mean_returns('FutRet') ## double check data integrity vs. external data
_ = pe.check_forward_and_mean_returns('DSpot')
_ = pe.check_forward_and_mean_returns('xomRet')
_ = pe.check_forward_and_mean_returns('bpRet')
_ = pe.check_forward_and_mean_returns('rdsaRet')
# %% port sims
ld,strat,srl,sr = pe.port_test('text','FutRet','FutRet',thresh=0,saveout=saveout)
ld,strat,srl,sr = pe.port_test('text','DSpot','FutRet',thresh=0,saveout=saveout)
ld,strat,srl,sr = pe.port_test('text','bpRet','bpRet',thresh=0,saveout=saveout)
ld,strat,srl,sr = pe.port_test('text','rdsaRet','rdsaRet',thresh=0,saveout=saveout)
ld,strat,srl,sr = pe.port_test('text','xomRet','xomRet',thresh=0,saveout=saveout)
# %% one-off cagr signal
_ = pe.port_test('text','FutRet','FutRet',signal_type='cagr_l1')
_ = pe.port_test('text','FutRet','FutRet',signal_type='cagr_l3')
_ = pe.port_test('text','FutRet','FutRet',signal_type='cagr_l12')
# %% compare text signal to momentum signals
pe.port_test_signals('text','FutRet','FutRet',saveout=saveout)
pe.port_test_signals('text','bpRet','bpRet',saveout=saveout)
pe.port_test_signals('text','rdsaRet','rdsaRet',saveout=saveout)
pe.port_test_signals('text','xomRet','xomRet',saveout=saveout)
# %% matrix
res = pe.port_test_matrix('text','pred',saveout=saveout)
res = pe.port_test_matrix('nontext','pred',saveout=saveout)
res = pe.port_test_matrix('both','pred',saveout=saveout)
res = pe.port_test_matrix('text','mean',saveout=saveout)
## momentum strategies
res = pe.port_test_matrix('text','cagr_l3',saveout=saveout)
res = pe.port_test_matrix('text','cagr_l6',saveout=saveout)
res = pe.port_test_matrix('text','cagr_l12',saveout=saveout)
res = pe.port_test_matrix('text','cagr_l12m1',saveout=saveout)
# %% double check some results from paper
## a Figure 9 plot
thed, oosr2s = pe.calc_selection_oosR2('text','FutRet','zero')
thed, oosr2s = pe.calc_selection_oosR2('text','xomRet','zero')
thed, oosr2s = pe.calc_selection_oosR2('text','rdsaRet','zero')
thed, oosr2s = pe.calc_selection_oosR2('text','DProd','zero')

## confirming the results in Table 6
eo = pe.selection_and_evaluation_oosR2('text','rolling')
eo = pe.selection_and_evaluation_oosR2('text','zero')
eo = pe.selection_and_evaluation_oosR2('both','zero')
# %% rdsa example
thed = pd.read_csv('c:/users/harry/code/energy/outofsample/res_Forward_lasso/0-1_rdsaRet.csv',
                   parse_dates=True,index_col=0)

## 2010-01-01 cutoff
used = thed[(thed.index >= '2010-01-01') & (thed['lookback_tercile_lag_4.0yr']==2)]
print(np.round(100 - 100*((used.true-used.pred)**2).sum()/((used.true-used['mean'])**2).sum(),2))

## 2009-11-06 cutoff
used = thed[(thed.index > '2009-11-06') & (thed['lookback_tercile_lag_4.0yr']==2)]
print(np.round(100 - 100*((used.true-used.pred)**2).sum()/((used.true-used['mean'])**2).sum(),2))