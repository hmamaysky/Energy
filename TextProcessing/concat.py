#!/user/kh3191/.conda/envs/nlp/bin/python

import numpy as np
import pandas as pd
from scipy import sparse

import glob
from tqdm import tqdm

# concat all info
folder_path = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info'
csv_files = glob.glob(folder_path + '/*_info.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)

output_file = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv'
concatenated_df.to_csv(output_file, index=False)

# concat dtm
folder_path = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C'
csv_files = glob.glob(folder_path + '/*_dtm.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)
output_file = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/dtm_concatenate.npz'
#concatenated_df.iloc[:,2:].to_csv(output_file, index=False)
print('Converting to sparse matrix')
sparse_matrix = sparse.csr_matrix(concatenated_df.iloc[:,2:])
print('Writing concatenated df')
np.savez(output_file, data=sparse_matrix.data, indices=sparse_matrix.indices,
         indptr=sparse_matrix.indptr, shape=sparse_matrix.shape)
