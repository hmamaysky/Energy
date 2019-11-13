#!/apps/anaconda2/bin/python

"""
    Program            : refer to run_info.sh 
    Function           : This code concatenate all the article measures and change the time to NY
"""

import pandas as pd
import sys
from dateutil import tz
from datetime import datetime


#################
####FUNCTIONS####
#################
def UTC_to_NY(row):
    x = row['TimeStamp']
    # Hardcode zones:
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/New_York')
    
    utc = datetime.strptime(x[0:19], '%Y-%m-%dT%H:%M:%S')
    
    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    utc = utc.replace(tzinfo=from_zone)
    
    # Convert time zone
    est = utc.astimezone(to_zone)
    est_timestamp = est.strftime('%Y-%m-%dT%H:%M:%S')+x[19:]
    return est_timestamp
################




aa = sys.argv[1]

path_sent = '/NOBACKUP/scratch/ra2826/oil-project/sentiment' 
path_topic = '/NOBACKUP/scratch/ra2826/oil-project/topic_allocation'
path_entropy = '/NOBACKUP/scratch/ra2826/oil-project/entropy'
path_total = '/NOBACKUP/scratch/ra2826/oil-project/total'
path_info = '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS'
path_out = '/NOBACKUP/scratch/ra2826/oil-project/info'                 


df0=pd.DataFrame()


df_sent = pd.read_csv(path_sent+ '/' + aa + 'sent.csv', delimiter=',')
df_sent.rename(columns={'sent': 'sentiment'}, inplace=True)
df_topic = pd.read_csv(path_topic+ '/' + aa + 'topic_alloc.csv', delimiter=',')
df_topic.drop(['headline'], axis=1, inplace=True)
df_entropy = pd.read_csv(path_entropy+ '/' + aa + '_entropy.csv', delimiter=',')
df_total = pd.read_csv(path_total+ '/' + aa + 'total.csv', delimiter=',')
df_info = pd.read_csv(path_info + '/' + 'oil_RTRS_' + aa + '.csv', delimiter=',')
df_info.drop(['augbod'], axis=1, inplace=True)


df_all0=df_info.join(df_sent['sentiment'])
df_all1=df_all0.join(df_topic)
df_all2=df_all1.join(df_entropy['entropy'])
df_all3=df_all2.join(df_total['total'])
df0=df_all3
        



df0['TimeStamp_NY'] = df0.apply(UTC_to_NY, axis=1)
df0.rename(columns={'TimeStamp': 'TimeStamp_UTC'}, inplace=True)

cols = ['Id', 'TimeStamp_UTC', 'TimeStamp_NY', 'subject', 'headline', 'entropy', 'total', 'sentiment', 'Topic1', 'Topic2', 'Topic3', 'Topic4', 'Topic5', 'Topic6', 'Topic7']
df0 = df0[cols]
        
df0.to_csv(path_out+ '/' + aa + '_info.csv',index=False)	





