#!/apps/anaconda2/bin/python

"""
    Program            : refer to run_topic.sh 
    Function           : This code calculates the allocation of the topics for each article 
"""

import string
import numpy as np
import pandas as pd
import sys
import os
import csv

j=sys.argv[1]
aa = j[-13:-7]

inputpath_dtm = '/NOBACKUP/scratch/ra2826/oil-project/dtm_Clustering_C'
inputpath_info = '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS'
FS='/user/user1/ra2826/oil_project/article_measures/topic_allocation'


df_dtm = pd.read_csv(inputpath_dtm + '/' + aa + 'dtm.csv', delimiter=',')

df_info = pd.read_csv(inputpath_info + '/' + 'oil_RTRS_' + aa + '.csv', delimiter=',')

topics = pd.read_csv(FS+'/'+'clustering_C.csv',sep=',',index_col=0)

topics_list=[]
for i in range(7):
    topics_list.append(topics.index[topics['Topic'] == i+1].tolist())


df0=pd.DataFrame()
for i in range(7):
    data=[]
    for index, row in df_dtm.iterrows():
        print(index)
        x=sum(df_dtm.loc[index,topics_list[i]])
        data.append(x)
    df0['Topic'+str(i+1)]=data


df0['sum'] = df0.sum(axis=1)
df0 = df0.loc[:,'Topic1':'Topic7'].div(df0['sum'], axis=0)
df0=df0.fillna(0)


df0['headline']=df_info['headline']


outputpath = '/NOBACKUP/scratch/ra2826/oil-project/topic_allocation'
df0.to_csv(outputpath + '/' + aa + 'topic_alloc.csv',index=False)		


