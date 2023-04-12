#!/user/kh3191/.conda/envs/nlp/bin/python

import numpy as np
import pandas as pd
import os

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--oldPath', type=str, 
       default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2020-11-16')
    parser.add_argument('--newPath', type=str, 
       default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023')
    parser.add_argument('--check', type=str, default='topic')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


if opt.check.lower()[:5] == 'topic':
    
    col_names = [f'Topic{i}' for i in range(1,8)]
    print('---- Correlation should be close to 1 if outputs are reliable.')
    # mismatch = []
    file_list = os.listdir(f'{opt.oldPath}/DataProcessing/topic_allocation')
    file_list.sort()
    for file in file_list:
        YYYYMM = file[:6]
        old = pd.read_csv(f'{opt.oldPath}/DataProcessing/topic_allocation/{file}')
        new = pd.read_csv(f'{opt.newPath}/DataProcessing/article_measure/topic_allocation/{YYYYMM}_topic_alloc.csv')

        old_stack = old[col_names].stack().reset_index(drop=True)
        new_stack = new[col_names].stack().reset_index(drop=True)
        try:
            corr = np.corrcoef(old_stack, new_stack).min()
            print(f'{YYYYMM}: {corr:.2f}')
        except ValueError:
            # mismatch.append(YYYYMM)
            print(f'{YYYYMM}: mismatched')