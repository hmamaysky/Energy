#!/user/kh3191/.conda/envs/nlp/bin/python

import pandas as pd
from scipy import sparse

import glob
from tqdm import tqdm


import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--combinedInfoPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info')
    parser.add_argument('--dtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)

# concat all info
csv_files = glob.glob(opt.combinedInfoPath + '/*_info.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)

output_file = f'{opt.outputPath}/info_concatenate.csv'
concatenated_df.to_csv(output_file, index=False)

# concat dtm
csv_files = glob.glob(opt.dtmPath + '/*_dtm.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)
output_file = f'{opt.outputPath}/dtm_concatenate.npz'
print('Converting to sparse matrix')
sparse_matrix = sparse.csr_matrix(concatenated_df.iloc[:,2:])
print('Writing concatenated df')
sparse.save_npz(output_file, sparse_matrix)
