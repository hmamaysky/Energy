#!/user/kh3191/.conda/envs/nlp/bin/python

"""
    Function           : This code combine all the article measures and change the time to NY
"""

import os
import pandas as pd
from dateutil import tz
from datetime import datetime
from tqdm import tqdm

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    parser.add_argument('--measurePath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info')
    parser.add_argument('--local_topic_model', type=bool, 
           default=False)
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


#################
####FUNCTIONS####
#################
# Hardcode zones:
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('America/New_York')
    
def UTC_to_NY(row):
    x = row['TimeStamp']
    utc = datetime.strptime(x[0:19], '%Y-%m-%dT%H:%M:%S')

    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    utc = utc.replace(tzinfo=from_zone)

    # Convert time zone
    est = utc.astimezone(to_zone)
    est_timestamp = est.strftime('%Y-%m-%dT%H:%M:%S')+x[19:]
    return est_timestamp
################


if __name__ == "__main__":

    if opt.local_topic_model:
        topic_alloc_files = os.listdir(f'{opt.measurePath}/rolling_topic_allocation')
        YYYYMM_start_list = [file[-29:-23] for file in topic_alloc_files]
    else:
        topic_alloc_files = os.listdir(f'{opt.measurePath}/topic_allocation')
    YYYYMM_list = [file[-22:-16] for file in topic_alloc_files]
        
    for k, YYYYMM_end in enumerate(tqdm(YYYYMM_list)):
        
        df_info = pd.read_csv(f"{opt.inputPath}/oil_{YYYYMM_end}_info.csv", delimiter=',')
        df_info.drop(['augbod'], axis=1, inplace=True)
        
        df_sent = pd.read_csv(f"{opt.measurePath}/sentiment/{YYYYMM_end}_sent.csv", delimiter=',')
        df_sent.rename(columns={'sent': 'sentiment'}, inplace=True)
        
        if opt.local_topic_model:
            YYYYMM_start = YYYYMM_start_list[k]
            df_topic = pd.read_csv(f"{opt.measurePath}/rolling_topic_allocation/{YYYYMM_start}_{YYYYMM_end}_topic_alloc.csv",
                                   delimiter=',')
        else:
            df_topic = pd.read_csv(f"{opt.measurePath}/topic_allocation/{YYYYMM_end}_topic_alloc.csv", delimiter=',')
        df_topic.drop(['headline'], axis=1, inplace=True)
        
        df_entropy = pd.read_csv(f"{opt.measurePath}/entropy/{YYYYMM_end}_entropy.csv", delimiter=',')
        df_total = pd.read_csv(f"{opt.measurePath}/total/{YYYYMM_end}_total.csv", delimiter=',')
    
        df = df_info.join(df_sent['sentiment'])\
                    .join(df_topic)\
                    .join(df_entropy['entropy'])\
                    .join(df_total['total'])
        
        df['TimeStamp_NY'] = df.apply(UTC_to_NY, axis=1)
        df.rename(columns={'TimeStamp': 'TimeStamp_UTC'}, inplace=True)

        cols = ['Id', 'TimeStamp_UTC', 'TimeStamp_NY', 'subject', 'headline', 'entropy', 'total', 'sentiment'] + \
                list(df_topic.columns)
        df = df[cols]

        if opt.local_topic_model:
            df.to_csv(f"{opt.outputPath}/{YYYYMM_start}_{YYYYMM_end}_info.csv",index=False)
        else:
            df.to_csv(f"{opt.outputPath}/{YYYYMM_end}_info.csv",index=False)
        
        
