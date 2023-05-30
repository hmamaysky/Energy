#!/user/kh3191/.conda/envs/nlp/bin/python
 
"""
    Function           : Input: info, 3gram, and 4gram files
                         Output: entropy files   
"""

import numpy as np
import pandas as pd
import torch
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=False)

from utils import get_clean4, NGRAM

import os
import glob

from tqdm import tqdm


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
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


if __name__ == "__main__":
    
    data = glob.glob(opt.inputPath + '/*.csv')
    data.sort()
    ngr3 = glob.glob(opt.ngPath + '/3gram/*.csv')
    ngr3.sort()
    ngr4 = glob.glob(opt.ngPath + '/4gram/*.csv')
    ngr4.sort()
    
    stop3_dicts = [get_stop_dict(ngr_file) for ngr_file in tqdm(ngr3)]
    torch.save(stop3_dicts, 'stop3_dicts.pt')
    stop4_dicts = [get_stop_dict(ngr_file) for ngr_file in tqdm(ngr4)]
    torch.save(stop4_dicts, 'stop4_dicts.pt')
    
    for fnum, file in tqdm(enumerate(data)):

        Temp = pd.read_csv(file,sep=',')

        # fnum 27 indicates the number of months to use as the trial data to use for entropy calculations
        if fnum >= 27:
            
            Temp['gram4'] = Temp['augbod'].parallel_apply(get_clean4)
            ngram4s = NGRAM()
            for f in Temp['gram4']:
                ngram4s.add_doc(f)
                      
#             import time
#             start_time = time.time()
            
            stop3_dict = get_stop_dict(ngr3[fnum-27:fnum-3])
            stop4_dict = get_stop_dict(ngr4[fnum-27:fnum-3])
#             print("--- %s seconds ---" % (time.time() - start_time))

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
