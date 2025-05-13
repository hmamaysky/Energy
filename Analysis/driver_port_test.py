from Energy.Analysis import port_test as pt
# %% read data
pe = pt.PortEngine()
pe.check_true_rets('FutRet')
pe.check_true_rets('xomRet')
pe.check_true_rets('bpRet')
pe.check_true_rets('rdsaRet')
# %% double check some results from paper
## a Figure 9 plot
pe.select_oosR2('text','FutRet')

## confirming the results in Table 6
eo = pe.eval_oosR2(type='text', lag=3, tercile=2, weight=1, var='xomRet')
# %% port sims
ld,wts = pe.port_test('text','DSpot','FutRet',lag=4,tercile=2,thresh=0)
ld,wts = pe.port_test('text','DSpot','bpRet',lag=4,tercile=2,thresh=0)
ld,wts = pe.port_test('text','DSpot','rdsaRet',lag=4,tercile=2,thresh=0)
ld,wts = pe.port_test('text','DSPot','xomRet',lag=4,tercile=2,thresh=0)
# %% matrix
forecast_vars = ['FutRet','bpRet','rdsaRet','xomRet','DSpot']
dep_vars = ['FutRet','bpRet','rdsaRet','xomRet']
res = pd.DataFrame(np.nan,index=forecast_vars,columns=dep_vars)
for forecast in forecast_vars:
    for dep in dep_vars:
        _, sr = pe.port_test('text',forecast,dep,4,2,0)
        res.loc[forecast,dep] = sr
        
print(res)
