#!/user/kh3191/.conda/envs/tmle/bin/python

import numpy as np
import pandas as pd
import itertools
import os
import torch

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--tagPath', type=str, 
                        default='energytag.csv')
    parser.add_argument('--oldPath', type=str, 
       default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2020-11-16')
    parser.add_argument('--newPath', type=str, 
       default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023')
    parser.add_argument('--check', type=str, default='topic')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


def flatten(lists, filterN2=True):
    
    set_tags = set(list(itertools.chain.from_iterable(lists)))
    if filterN2:
        return {tag for tag in set_tags if tag.startswith('N2')}
    else:
        return set_tags


if opt.check.lower()[:5] == 'topic':
    
    col_names = [f'Topic{i}' for i in range(1,8)]
    print('---- Correlation should be close to 1 if outputs are reliable.')

    file_list = os.listdir(f'{opt.oldPath}/DataProcessing/topic_allocation')
    file_list.sort()
    
    tags = pd.read_csv(opt.tagPath)['energytag'].values
    missingTags = {}
    
    for file in file_list:
        YYYYMM = file[:6]
        old = pd.read_csv(f'{opt.oldPath}/DataProcessing/topic_allocation/{file}')
        new = pd.read_csv(f'{opt.newPath}/DataProcessing/article_measure/topic_allocation/{YYYYMM}_topic_alloc.csv')

        if len(old) == len(new):
            old_stack = old[col_names].stack().reset_index(drop=True)
            new_stack = new[col_names].stack().reset_index(drop=True)
            corr = np.corrcoef(old_stack, new_stack).min()
            print(f'{YYYYMM}: {corr:.3f}')
            
            
        else:
            raw = pd.read_csv(f'{opt.newPath}/DataProcessing/info/{YYYYMM}_info.csv')
            raw = raw.sort_values('TimeStamp').groupby('PNAC').first().reset_index()
            
            old_selected = raw['headline'].isin(old['headline'])
            list_tags = raw.loc[old_selected, 'subject'].apply(eval).to_list()
            list_tags = flatten(list_tags)
            not_list_tags = raw.loc[~old_selected, 'subject'].apply(eval).to_list()
            not_list_tags = flatten(not_list_tags)
            
            diff_list_tags = list_tags.difference(not_list_tags)
            diff_list_tags = {i[3:] for i in diff_list_tags}
            diff_list_tags = diff_list_tags.difference(tags)
            print(f'{YYYYMM} missing tags: {len(diff_list_tags)}')
            missingTags[YYYYMM] = diff_list_tags
            
    torch.save(missingTags, 'missingTags.pt')
            