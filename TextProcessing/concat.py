#!/user/kh3191/.conda/envs/nlp/bin/python

import glob
import pandas as pd
from tqdm import tqdm

# concat all info
folder_path = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info'
csv_files = glob.glob(folder_path + '/*.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)

output_file = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv'
concatenated_df.to_csv(output_file, index=False)

# concat dtm
folder_path = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C'
csv_files = glob.glob(folder_path + '/*.csv')
csv_files.sort()

concatenated_df = pd.concat([pd.read_csv(file) for file in tqdm(csv_files)], ignore_index=True)
output_file = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C/dtm_concatenate.csv'
print('Writing concatenated df')
concatenated_df.iloc[:,2:].to_csv(output_file, index=False)
