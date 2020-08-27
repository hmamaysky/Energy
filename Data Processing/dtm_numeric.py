#!/apps/anaconda2/bin/python
 
"""
    Program            : refer to run_dtm_numeric.sh 
    Function           : This code converts the words and article ids in dtm files to numbers.
    					The reason is: this is the format that the sparse matrix in cosine code gets	
"""

import pandas as pd
import sys


j= sys.argv[1]
dtmpath= '/work/hw2676/Energy/dtm_Clustering_C'
df = pd.read_csv(dtmpath + '/' + j)
df = pd.melt(df, id_vars=list(df)[:2], value_vars=list(df)[2:], var_name='words', value_name='freq')
df = df[df['freq']!=0]  

word_column = df.words.tolist()
id_column = df.Id.tolist()
freq_column = df.freq.tolist()


words_test = pd.read_csv('/user/hw2676/code/Energy/DataProcessing/article_measure/dtm/clustering_C.csv', sep=',',names=['word','Topic','freq'])
word_set = words_test['word'][1:].tolist()


idpath = '/work/hw2676/Energy/concat'
df_selected = pd.read_csv(idpath + '/' + 'Info_concat.csv',sep=',')
id_set = df_selected.Id.tolist()


dic_word = {k: v for v, k in enumerate(word_set)}
def func_word(item):
    return dic_word[item]  


dic_id = {k: v for v, k in enumerate(id_set)}
def func_id(item):
    return dic_id[item]  


V = freq_column
I = map(func_id,id_column)
J = map(func_word,word_column)


df_out = pd.DataFrame()
df_out['freq'] = V
df_out['Id'] = I
df_out['words'] = J

outputpath= '/work/hw2676/Energy/dtm_numeric'
df_out.to_csv(outputpath + '/' + j,index=False)


