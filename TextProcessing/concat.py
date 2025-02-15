#!/user/kh3191/.conda/envs/nlp/bin/python

import pandas as pd
from scipy import sparse

from glob import glob
from tqdm import tqdm
# ./concat.py --concat_info='' --concat_dtm='' --dtmPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_Clustering_C --save_monthly_freq=True

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--combinedInfoPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info')
    parser.add_argument('--concat_info', type=bool, 
           default=True)
    parser.add_argument('--concat_dtm', type=bool, 
           default=True)
    parser.add_argument('--dtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat')
    parser.add_argument('--save_monthly_freq', type=bool, 
           default=False)
    parser.add_argument('--freqPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_freq')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)

if opt.concat_info:
    # concat all info
    info_files = glob(opt.combinedInfoPath + '/*_info.csv')
    info_files.sort()

    concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(info_files)], ignore_index=True)
    output_file = f'{opt.outputPath}/info_concatenate.csv'
    concatenated_df.to_csv(output_file, index=False)


if opt.concat_dtm:
    # concat dtm
    dtm_files = glob(opt.dtmPath + '/*_dtm.csv')
    dtm_files.sort()
    
    concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(dtm_files)], ignore_index=True)
    output_file = f'{opt.outputPath}/dtm_concatenate.npz'
    print('Converting to sparse matrix')
    sparse_matrix = sparse.csr_matrix(concatenated_df.iloc[:,2:])
    print('Writing concatenated df')
    sparse.save_npz(output_file, sparse_matrix)

    
if opt.save_monthly_freq:
    # save word freq at monthly level
    dtm_files = glob(opt.dtmPath + '/*_dtm.csv')
    dtm_files.sort()
    
    for file in tqdm(dtm_files):
        YYYYMM = file[-14:-8]
        freq_series = pd.read_csv(file).iloc[:,2:].sum()
        freq_series.to_csv(f'{opt.freqPath}/{YYYYMM}_freq.csv')