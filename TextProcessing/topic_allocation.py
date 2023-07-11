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
       default='clustering_C.csv')
    parser.add_argument('--inputPath_info', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    parser.add_argument('--inputPath_dtm', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C')
    parser.add_argument('--outputPath', type=str, 
       default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/topic_allocation')
    parser.add_argument('--n_topics', type=int, default=7)
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


for file in tqdm(os.listdir(opt.inputPath_dtm)):
    YYYYMM = file[-15:-8]

    df_dtm = pd.read_csv(f'{opt.inputPath_dtm}/{YYYYMM}_dtm.csv', delimiter=',')
    df_info = pd.read_csv(f'{opt.inputPath_info}/oil_{YYYYMM}_info.csv', delimiter=',')
    df_topics = pd.read_csv(opt.inputWordsPath, sep=',', index_col=0)
    topics_list = [df_topics.index[df_topics['Topic'] == i].tolist() for i in range(1,opt.n_topics+1)]

    df0 = pd.DataFrame()
    for i, topics in enumerate(topics_list):
        df0['Topic'+str(i+1)] = df_dtm[topics].sum(axis=1)

    df0['sum'] = df0.sum(axis=1)
    df0 = df0.loc[:,'Topic1':f'Topic{opt.n_topics}'].div(df0['sum'], axis=0)
    df0 = df0.fillna(0)

    df0['headline']=df_info['headline']
    df0.to_csv(f'{opt.outputPath}/{YYYYMM}_topic_alloc.csv',index=False)

