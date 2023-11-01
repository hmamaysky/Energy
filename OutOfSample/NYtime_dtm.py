#!/user/kh3191/.conda/envs/nlp/bin/python
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

from glob import glob
from tqdm import tqdm
#tqdm.pandas()
from info import UTC_to_NY
from date_fixed_measures import oil_date

if __name__ == '__main__':
    
    dtmPath = '/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C'
    dtm_files = glob(f'{dtmPath}/*.csv')
    dtm_files.sort()
    
    df = pd.concat([pd.read_csv(file) for file in tqdm(dtm_files)])
    df.loc[:,'TimeStamp'] = df.parallel_apply(UTC_to_NY, axis=1)
    df['TimeStamp_date'] = df['TimeStamp'].parallel_apply(lambda sample: oil_date(sample, store_date=True))
    df['TimeStamp_date'] = pd.to_datetime(df['TimeStamp_date'])
    # groupby NYtime year and month
    groups = df.groupby([df['TimeStamp_date'].dt.year, df['TimeStamp_date'].dt.month])

    for (year, month), group in tqdm(groups):
        group.drop(columns=['TimeStamp_date'], inplace=True)
        filename = f"/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_Clustering_C/{year}{month:02}_dtm.csv"
        group.to_csv(filename, index=False)
