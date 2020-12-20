#!/apps/anaconda2/bin/python


#This code gets the daily aggregates 
import pandas as pd

##########################
#Read the input

def Unclassified_func(s):
    if s>.98: return 0
    else: return 1
    

wkdir = '/work/hw2676/Energy/concat/'

df = pd.read_csv(wkdir + 'date_fixed_article_level_measures.csv',sep=',')
df['sum'] = df['Topic1'] + df['Topic2'] + df['Topic3'] + df['Topic4'] + df['Topic5'] + df['Topic6'] + df['Topic7']
df['Unclassified'] = df['sum'].apply(Unclassified_func)



df=df.dropna()
##########################

#one measure count of articles
df_count=df.groupby('date').size()

#one entropy measure
df_ent=df.groupby('date').apply(lambda dfx: (dfx['entropy'] * dfx['total']).sum() / dfx['total'].sum())

#six topic measures
df_t1=df.groupby('date').apply(lambda dfx: (dfx['Topic1'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t2=df.groupby('date').apply(lambda dfx: (dfx['Topic2'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t3=df.groupby('date').apply(lambda dfx: (dfx['Topic3'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t4=df.groupby('date').apply(lambda dfx: (dfx['Topic4'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t5=df.groupby('date').apply(lambda dfx: (dfx['Topic5'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t6=df.groupby('date').apply(lambda dfx: (dfx['Topic6'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t7=df.groupby('date').apply(lambda dfx: (dfx['Topic7'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t8=df.groupby('date').apply(lambda dfx: (dfx['Unclassified'] * dfx['total']).sum() / float(dfx['total'].sum()))


#six topic-sentiment measures
df_t1_s=df.groupby('date').apply(lambda dfx: (dfx['Topic1'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t2_s=df.groupby('date').apply(lambda dfx: (dfx['Topic2'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t3_s=df.groupby('date').apply(lambda dfx: (dfx['Topic3'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t4_s=df.groupby('date').apply(lambda dfx: (dfx['Topic4'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t5_s=df.groupby('date').apply(lambda dfx: (dfx['Topic5'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t6_s=df.groupby('date').apply(lambda dfx: (dfx['Topic6'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t7_s=df.groupby('date').apply(lambda dfx: (dfx['Topic7'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))

df_t8_s=df.groupby('date').apply(lambda dfx: (dfx['Unclassified'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))



d = {'article count' : df_count, 'entropy': df_ent, 
     'Topic 1' : df_t1, 'Topic 2' : df_t2,'Topic 3' : df_t3, 
     'Topic 4' : df_t4, 'Topic 5' : df_t5, 'Topic 6' : df_t6, 'Topic 7' : df_t7,
     'Unclassified' : df_t8,
     'Topic-Sentiment 1' : df_t1_s, 'Topic-Sentiment 2' : df_t2_s,'Topic-Sentiment 3' : df_t3_s, 
     'Topic-Sentiment 4' : df_t4_s, 'Topic-Sentiment 5' : df_t5_s, 'Topic-Sentiment 6' : df_t6_s, 'Topic-Sentiment 7' : df_t7_s,
     'Unclassified-Sentiment' : df_t8_s}
df_out = pd.DataFrame(d)

cols = ['article count','entropy',
     'Topic 1','Topic 2','Topic 3', 
     'Topic 4','Topic 5','Topic 6','Topic 7','Unclassified',
     'Topic-Sentiment 1','Topic-Sentiment 2','Topic-Sentiment 3',
     'Topic-Sentiment 4','Topic-Sentiment 5','Topic-Sentiment 6','Topic-Sentiment 7',
     'Unclassified-Sentiment']
df_out = df_out[cols]
df_out.to_csv('NYtime_daily_level_measures_C_2020-05-03.csv')