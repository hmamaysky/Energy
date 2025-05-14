from Energy.Analysis import port_test as pt
# %% read data
pe = pt.PortEngine()
pe.check_forward_and_mean_returns('FutRet')
pe.check_forward_and_mean_returns('xomRet')
pe.check_forward_and_mean_returns('bpRet')
pe.check_forward_and_mean_returns('rdsaRet')
# %% double check some results from paper
## a Figure 9 plot
pe.select_oosR2('text','FutRet')

## confirming the results in Table 6
eo = pe.eval_oosR2(type='text', lag=3, tercile=2, weight=1, var='xomRet')
# %% port sims
ld,wts,sr = pe.port_test('text','FutRet','FutRet',lag=4,tercile=2,thresh=0)
ld,wts,sr = pe.port_test('text','bpRet','bpRet',lag=4,tercile=2,thresh=0)
ld,wts,sr = pe.port_test('text','rdsaRet','rdsaRet',lag=4,tercile=2,thresh=0)
ld,wts,sr = pe.port_test('text','xomRet','xomRet',lag=4,tercile=2,thresh=0)
# %% matrix
pe.port_test_matrix('pred')