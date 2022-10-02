from Energy.Analysis import oos_energy as oos
pd.set_option('display.max_columns', 150)
# %% get the data
oo = oos.OOSAnalysis()
# %% replicate results in paper
oo.gen_old_table('1_1')
oo.gen_old_table('2_2')
# %% do predictions versus
oo.blended_oos('1_1','base')
oo.blended_oos('1_1','text')
oo.blended_oos('1_1','full')
oo.blended_oos('2_2','base')
oo.blended_oos('2_2','text')
oo.blended_oos('2_2','full')