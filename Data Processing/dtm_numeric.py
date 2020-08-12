#!/apps/anaconda2/bin/python
 
"""
    Program            : refer to run_dtm_numeric.sh 
    Function           : This code converts the words and article ids in dtm files to numbers.
    					The reason is: this is the format that the sparse matrix in cosine code gets	
"""

import pandas as pd
import sys


j= sys.argv[1]
dtmpath= '/NOBACKUP/scratch/ra2826/oil-project/dtm_Clustering_C'
df = pd.read_csv(dtmpath + '/' + j)
df = pd.melt(df, id_vars=list(df)[:2], value_vars=list(df)[2:], var_name='words', value_name='freq')
df = df[df['freq']!=0]

word_column = df.words.tolist()
id_column = df.Id.tolist()
freq_column = df.freq.tolist()


words_test = pd.read_csv('/user/user1/ra2826/oil_project/article_measures/dtm/clustering_C.csv', sep=',')
word_set = words_test.word.tolist()


idpath = '/NOBACKUP/scratch/ra2826/oil-project/concat'
df_selected = pd.read_csv(idpath + '/' + 'oil_articles.csv',sep=',')
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

outputpath= '/NOBACKUP/scratch/ra2826/oil-project/dtm_numeric'
df_out.to_csv(outputpath + '/' + j,index=False)


