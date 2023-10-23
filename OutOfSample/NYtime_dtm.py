#!/user/kh3191/.conda/envs/nlp/bin/python
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=False)

from glob import glob
from tqdm import tqdm
#tqdm.pandas()
from info import UTC_to_NY

if __name__ == '__main__':
    
    dtmPath = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C'
    dtm_files = glob(f'{dtmPath}/*.csv')
    dtm_files.sort()
    
    df = pd.concat([pd.read_csv(file) for file in tqdm(dtm_files)])
    df.loc[:,'TimeStamp'] = df.parallel_apply(lambda row: UTC_to_NY(row, to_str=False), axis=1)
    df['TimeStamp'] = pd.to_datetime(df['TimeStamp'])
    # groupby NYtime year and month
    groups = df.groupby([df['TimeStamp'].dt.year, df['TimeStamp'].dt.month])

    for (year, month), group in tqdm(groups):
        filename = f"/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_Clustering_C/{year}{month:02}_dtm.csv"
        group.to_csv(filename, index=False)
