#!/user/kh3191/.conda/envs/nlp/bin/python

"""
    Function           : This code calculates the allocation of the topics for each article 
"""

import numpy as np
import pandas as pd
import os
from tqdm import tqdm

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputWordsPath', type=str, 
       default='clustering_C.csv',
       help='Clustering csvs are used to determine which words are contained in each topic.')
    parser.add_argument('--inputPath_info', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info',
       help='Info files are used to retrieve headlines.')
    parser.add_argument('--inputPath_dtm', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C',
       help='Dtms are used to sum up total word frequencies in each topic.')
    parser.add_argument('--outputPath', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/topic_allocation')
    parser.add_argument('--local_topic_model', type=bool, 
       default=False)
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


if opt.local_topic_model:
    assert not opt.inputWordsPath.endswith('csv')
    inputWordsFiles = os.listdir(opt.inputWordsPath)
    inputWordsFiles.sort()
    YYYYMM_list = [file[-10:-4] for file in inputWordsFiles]
    YYYYMM_start_list = [file[-17:-11] for file in inputWordsFiles]
else:
    assert opt.inputWordsPath.endswith('csv')
    YYYYMM_list = [file[-14:-8] for file in os.listdir(opt.inputPath_dtm)]
    

for k, YYYYMM_end in enumerate(tqdm(YYYYMM_list)):
    
    if opt.local_topic_model:
        from concurrent.futures import ThreadPoolExecutor
        date_range_list = YYYYMM_start_list[k:YYYYMM_start_list.index(YYYYMM_end)]
        with ThreadPoolExecutor() as executor:
            dtm_frames = list(executor.map(lambda YYYYMM: pd.read_csv(f'{opt.inputPath_dtm}/{YYYYMM}_dtm.csv', 
                                                                      delimiter=','),
                                           date_range_list))
            df_dtm = pd.concat(dtm_frames)
        with ThreadPoolExecutor() as executor:
            info_frames = list(executor.map(lambda YYYYMM: pd.read_csv(f'{opt.inputPath_info}/oil_{YYYYMM}_info.csv',
                                                                       delimiter=','), 
                                            date_range_list))
            df_info = pd.concat(info_frames)
#         df_dtm = pd.concat([pd.read_csv(f'{opt.inputPath_dtm}/{YYYYMM}_dtm.csv', delimiter=',') 
#                             for YYYYMM in date_range_list])
#         df_info = pd.concat([pd.read_csv(f'{opt.inputPath_info}/oil_{YYYYMM}_info.csv', delimiter=',') 
#                              for YYYYMM in date_range_list])
        YYYYMM_start = YYYYMM_start_list[k]
        df_topics = pd.read_csv(f'{opt.inputWordsPath}/clustering_C_{YYYYMM_start}_{YYYYMM_end}.csv', sep=',', index_col=0)
    else:
        df_dtm = pd.read_csv(f'{opt.inputPath_dtm}/{YYYYMM_end}_dtm.csv', delimiter=',')
        df_info = pd.read_csv(f'{opt.inputPath_info}/oil_{YYYYMM_end}_info.csv', delimiter=',')
        df_topics = pd.read_csv(opt.inputWordsPath, sep=',', index_col=0)
        
    n_topics = df_topics['Topic'].max()
    topics_list = [df_topics.index[df_topics['Topic'] == i].tolist() for i in range(1,n_topics+1)]

    df0 = pd.DataFrame()
    for i, topics in enumerate(topics_list):
        df0['Topic'+str(i+1)] = df_dtm[topics].sum(axis=1)

    df0['sum'] = df0.sum(axis=1)
    df0 = df0.loc[:,'Topic1':f'Topic{n_topics}'].div(df0['sum'], axis=0)
    df0 = df0.fillna(0)

    df0['headline'] = df_info['headline']
    if opt.local_topic_model:
        df0.to_csv(f'{opt.outputPath}/{YYYYMM_start}_{YYYYMM_end}_topic_alloc.csv', index=False)
    else:
        df0.to_csv(f'{opt.outputPath}/{YYYYMM_end}_topic_alloc.csv', index=False)
