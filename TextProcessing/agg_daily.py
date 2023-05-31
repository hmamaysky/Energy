#!/user/kh3191/.conda/envs/nlp/bin/python


#This code gets the daily aggregates 
import pandas as pd
from tqdm import tqdm

##########################
#Read the input

wkdir = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat'
df = pd.read_csv(f"{wkdir}/date_fixed_article_level_measures.csv", sep=',')

topic_cols = [i.startswith('Topic') for i in df.columns]
n_topics = sum(topic_cols)
df['sum'] = df.loc[:,topic_cols].sum(axis=1)
df['Unclassified'] = (df['sum'] <= 0.98).astype(int)

df = df.dropna()
df_daily = df.groupby('date')
##########################

#one measure count of articles
df_count = df_daily.size()

#one entropy measure
df_ent = df_daily.apply(lambda dfx: (dfx['entropy'] * dfx['total']).sum() / dfx['total'].sum())

#topic measures
df_t1, df_t2, df_t3, df_t4, df_t5, df_t6, df_t7 = [df_daily.apply(
    lambda dfx: (dfx[f'Topic{i+1}'] * dfx['total']).sum() / float(dfx['total'].sum())
) for i in tqdm(range(7))]

df_t8 = df_daily.apply(lambda dfx: (dfx['Unclassified'] * dfx['total']).sum() / float(dfx['total'].sum()))


#topic-sentiment measures
df_t1_s, df_t2_s, df_t3_s, df_t4_s, df_t5_s, df_t6_s, df_t7_s = [df_daily.apply(
    lambda dfx: (dfx[f'Topic{i+1}'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum())
) for i in tqdm(range(7))]

df_t8_s = df_daily.apply(lambda dfx: (dfx['Unclassified'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))



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
df_out.to_csv('NYtime_daily_level_measures_C_2023.csv')
