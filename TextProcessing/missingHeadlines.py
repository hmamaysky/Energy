#!/user/kh3191/.conda/envs/tmle/bin/python

import pandas as pd
import os
from tqdm import tqdm

oldPath = '/shared/share_mamaysky-glasserman/energy_drivers/2020-11-16'
newPath = '/shared/share_mamaysky-glasserman/energy_drivers/2023'

file_list = os.listdir(f'{oldPath}/DataProcessing/topic_allocation')
file_list.sort()

save_dir = 'missingHeadlines.txt'

with open(save_dir, 'w') as f:
    f.write('')
        
for file in tqdm(file_list):
    YYYYMM = file[:6]
    old = pd.read_csv(f'{oldPath}/DataProcessing/topic_allocation/{file}')
    new = pd.read_csv(f'{newPath}/DataProcessing/article_measure/topic_allocation/{YYYYMM}_topic_alloc.csv')
    if len(old) != len(new):
        raw = pd.read_csv(f'{newPath}/DataProcessing/info/{YYYYMM}_info.csv')
        raw = raw.sort_values('TimeStamp').groupby('PNAC').first().reset_index()

        old_selected = raw['headline'].isin(old['headline'])
        new_selected = raw['headline'].isin(new['headline'])
        headlineDiff = raw.loc[(old_selected) & (~new_selected), 'headline'].to_list()
        with open(save_dir, 'a') as f:
            f.write('\n'.join(headlineDiff))
            f.write('\n')
