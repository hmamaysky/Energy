#!/user/kh3191/.conda/envs/nlp/bin/python
 
"""
    Function           : Input: info, 3gram, and 4gram files
                         Output: entropy files   
"""

import numpy as np
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=False)

from utils import get_clean4, NGRAM

import os
import glob

from tqdm import tqdm


######################################################## 
# 
# File Paths 
# 
########################################################
import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--ngPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure')
    parser.add_argument('--inputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    parser.add_argument('--monthTrials', type=int, default=27)
    parser.add_argument('--monthWindow', type=int, default=24)
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


# from collections import defaultdict
# def get_stop_dict(ngr_files):
#     stop_dict = defaultdict(int)
#     dfs = [pd.read_csv(l, index_col=0) for l in ngr_files]
#     df = pd.concat(dfs)
#     for word, freq in zip(df['word'], df['freq']):
#         stop_dict[word] += freq
#     return stop_dict

from collections import Counter
def get_stop_dict(ngr_file):
    df = pd.read_csv(ngr_file, index_col=0)
    stop_dict = Counter(dict(zip(df['word'], df['freq'])))
    return stop_dict

import json
def save_list_counters(counters, file_name):
    counters_data = [{k: v for k, v in counter.items()} for counter in counters]
    with open(file_name, 'w') as f:
        json.dump(counters_data, f)
        
def load_list_counters(file_name):
    with open(file_name, 'r') as f:
        counters_data = json.load(f)
        counters = [Counter(counter_dict) for counter_dict in counters_data]
    return counters

def get_rolling_counters(counters, fnum):
    stop_dict = Counter()
    for i in range(fnum-opt.monthTrials, fnum-opt.monthTrials+opt.monthWindow):
        stop_dict += counters[i]
    return stop_dict


if __name__ == "__main__":
    
    data = glob.glob(opt.inputPath + '/*.csv')
    data.sort()

    stop34_dicts = []
    for n in [3,4]:
        file_name = f"stop{n}_dicts.json"
        ngr = glob.glob(f"{opt.ngPath}/{n}gram/*.csv")
        ngr.sort()
        if not os.path.isfile(file_name):
            stop_dicts = [get_stop_dict(ngr_file) for ngr_file in tqdm(ngr)]
            save_list_counters(stop_dicts, file_name)
            stop34_dicts.append(stop_dicts)
        else:
            print(f"Loading {file_name}")
            stop34_dicts.append(load_list_counters(file_name))
    stop3_dicts, stop4_dicts = stop34_dicts
    
    for fnum, file in tqdm(enumerate(data)):

        Temp = pd.read_csv(file,sep=',')

        # fnum monthTrials indicates the number of months to use as the trial data to use for entropy calculations
        if fnum >= opt.monthTrials:
                      
            stop3_dict = get_rolling_counters(stop3_dicts, fnum)
            stop4_dict = get_rolling_counters(stop4_dicts, fnum)
            
            Temp['gram4'] = Temp['augbod'].parallel_apply(get_clean4)
            ngram4s = NGRAM()
            for f in Temp['gram4']:
                ngram4s.add_doc(f)
                
            entropy = []
            for a in ngram4s.sparse:
                a_sum = sum(a.values())
                p = [l/a_sum for l in a.values()]
                M = [(stop4_dict.get(l,0) + 1)/(stop3_dict.get('.'.join(l.split('.')[:3]),0) + 10) for l in a.keys()]
                m = [-pi * np.log(Mi) for pi, Mi in zip(p,M)]
                entropy.append(sum(m))
                
        else:
            entropy = [np.nan] * len(Temp)

        entpd = pd.DataFrame({'entropy':entropy})
        entpd = pd.concat([entpd,Temp[['Id']]],axis=1)
        
        YYYYMM = data[fnum][-15:-9]
        entpd.to_csv(f"{opt.ngPath}/entropy/{YYYYMM}_entropy.csv",index=False)
