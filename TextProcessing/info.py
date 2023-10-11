#!/user/kh3191/.conda/envs/nlp/bin/python

"""
    Function           : This code combine all the article measures and change the time to NY
"""

import os
import pandas as pd
from dateutil import tz
from datetime import datetime
from tqdm import tqdm
from utils import generate_month_list

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
    parser.add_argument('--local_topic_model', type=bool, default=False)
    parser.add_argument('--rolling_index', type=int, default=0)
    opt = parser.parse_args()
    return opt


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

def read_csv_file(opt, file_type, YYYYMM):
    
    if file_type == "info":
        path = f"{opt.inputPath}/oil_{YYYYMM}_info.csv"
        return pd.read_csv(path, delimiter=',').drop(['augbod'], axis=1)
    
    elif file_type == "sent":
        path = f"{opt.measurePath}/sentiment/{YYYYMM}_sent.csv"
        return pd.read_csv(path, delimiter=',').rename(columns={'sent': 'sentiment'})
    
    elif file_type == "topic":
        path = f"{opt.measurePath}/topic_allocation/{YYYYMM}_topic_alloc.csv"
        return pd.read_csv(path, delimiter=',').drop(['headline'], axis=1)
    
    elif file_type == "entropy":
        path = f"{opt.measurePath}/entropy/{YYYYMM}_entropy.csv"
        return pd.read_csv(path, delimiter=',')
    
    elif file_type == "total":
        path = f"{opt.measurePath}/total/{YYYYMM}_total.csv"
        return pd.read_csv(path, delimiter=',')
    
    else:
        raise ValueError("Invalid file type")

################


if __name__ == "__main__":
    
    opt = parse_option()
    print(opt)

    if opt.local_topic_model:
        topic_alloc_files = os.listdir(f'{opt.measurePath}/rolling_topic_allocation')
    else:
        topic_alloc_files = os.listdir(f'{opt.measurePath}/topic_allocation')
    topic_alloc_files.sort()
    YYYYMM_start_list = [file[-29:-23] for file in topic_alloc_files]
    YYYYMM_list = [file[-22:-16] for file in topic_alloc_files]
        
    file_types = ['info', 'sent', 'topic', 'entropy', 'total']
    
    for k in [opt.rolling_index]:
        assert 0 <= k <= 266
        YYYYMM_end = YYYYMM_list[k]
        
        if opt.local_topic_model:
            YYYYMM_start = YYYYMM_start_list[k]
            date_range_list = generate_month_list(YYYYMM_start, YYYYMM_end)
        
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor() as executor:
                info_frames = list(executor.map(lambda YYYYMM: read_csv_file(opt, "info", YYYYMM), date_range_list))
                sent_frames = list(executor.map(lambda YYYYMM: read_csv_file(opt, "sent", YYYYMM), date_range_list))
                entropy_frames = list(executor.map(lambda YYYYMM: read_csv_file(opt, "entropy", YYYYMM), date_range_list))
                total_frames = list(executor.map(lambda YYYYMM: read_csv_file(opt, "total", YYYYMM), date_range_list))

            df_info = pd.concat(info_frames).reset_index(drop=True)
            df_sent = pd.concat(sent_frames).reset_index(drop=True)
            df_topic = pd.read_csv(f"{opt.measurePath}/rolling_topic_allocation/{YYYYMM_start}_{YYYYMM_end}_topic_alloc.csv",
                                   delimiter=',').reset_index(drop=True)
            df_entropy = pd.concat(entropy_frames).reset_index(drop=True)
            df_total = pd.concat(total_frames).reset_index(drop=True)
        
        else:
            df_info, df_sent, df_topic, df_entropy, df_total = [read_csv_file(opt, 
                                                                              file_type, 
                                                                              YYYYMM_end) for file_type in file_types]
        
        df = pd.concat([df_info, df_sent['sentiment'], df_topic.iloc[:,:-1], df_entropy['entropy'], df_total['total']], axis=1)
        
        df['TimeStamp_NY'] = df.apply(UTC_to_NY, axis=1)
        df.rename(columns={'TimeStamp': 'TimeStamp_UTC'}, inplace=True)

        if opt.local_topic_model:
            cols = ['Id', 'TimeStamp_UTC', 'TimeStamp_NY', 'entropy', 'total', 'sentiment'] + list(df_topic.columns)
        else:
            cols = ['Id', 'TimeStamp_UTC', 'TimeStamp_NY', 'subject', 'headline', 'entropy', 'total', 'sentiment'] + \
                list(df_topic.columns)
        df = df[cols]

        if opt.local_topic_model:
            df.to_csv(f"{opt.outputPath}/{YYYYMM_start}_{YYYYMM_end}_info.csv",index=False)
        else:
            df.to_csv(f"{opt.outputPath}/{YYYYMM_end}_info.csv",index=False)
        
        
