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
df_t_list = [df_daily.apply(
    lambda dfx: (dfx[f'Topic{i+1}'] * dfx['total']).sum() / float(dfx['total'].sum())
) for i in tqdm(range(n_topics))]

df_t_unclass = df_daily.apply(lambda dfx: (dfx['Unclassified'] * dfx['total']).sum() / float(dfx['total'].sum()))


#topic-sentiment measures
df_t_s_list = [df_daily.apply(
    lambda dfx: (dfx[f'Topic{i+1}'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum())
) for i in tqdm(range(n_topics))]

df_t_s_unclass = df_daily.apply(lambda dfx: (dfx['Unclassified'] * dfx['sentiment'] * dfx['total']).sum() / float(dfx['total'].sum()))


d = {'article count' : df_count, 'entropy': df_ent, 
     'Unclassified' : df_t_unclass, 'Unclassified-Sentiment' : df_t_s_unclass}
for i, df_t in enumerate(df_t_list):
    d[f"Topic {i+1}"] = df_t
for i, df_t_s in enumerate(df_t_s_list):
    d[f"Topic-Sentiment {i+1}"] = df_t_s
df_out = pd.DataFrame(d)

# cols = ['article count','entropy'] + [f"Topic {i+1}" for i in range(n_topics)] + ['Unclassified'] +\
#         [f"Topic-Sentiment {i+1}" for i in range(n_topics)] + ['Unclassified-Sentiment']
# df_out = df_out[cols]
df_out.to_csv('NYtime_daily_level_measures_C_2023.csv')
