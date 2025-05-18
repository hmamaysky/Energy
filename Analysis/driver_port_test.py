from Energy.Analysis import port_test as pt
saveout = False
# %% read data
pe = pt.PortEngine()
pe.check_forward_and_mean_returns('FutRet')
pe.check_forward_and_mean_returns('xomRet')
pe.check_forward_and_mean_returns('bpRet')
pe.check_forward_and_mean_returns('rdsaRet')
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
# %% port sims
ld,wts,srl,sr = pe.port_test('text','FutRet','FutRet',thresh=0,saveout=saveout)
ld,wts,srl,sr = pe.port_test('text','bpRet','bpRet',thresh=0,saveout=saveout)
ld,wts,srl,sr = pe.port_test('text','rdsaRet','rdsaRet',thresh=0,saveout=saveout)
ld,wts,srl,sr = pe.port_test('text','xomRet','xomRet',thresh=0,saveout=saveout)
# %% matrix
res = pe.port_test_matrix('pred',saveout=saveout)
res = pe.port_test_matrix('mean',saveout=saveout)
