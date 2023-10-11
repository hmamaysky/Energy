#!/user/kh3191/.conda/envs/nlp/bin/python
"""
Adapted from Roya's codes
    Function           : This code fixes the dates on info files based on the oil price eastern closing time 
"""
import pandas as pd
import calendar
import datetime

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--concatPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat')
    parser.add_argument('--local_topic_model', type=bool, default=False)
    parser.add_argument('--rolling_index', type=int, default=0)
    opt = parser.parse_args()
    return opt

####################
##### FUNCTION #####
####################
def oil_date(sample):
    my_datetime = datetime.datetime.strptime(sample[0:19], '%Y-%m-%dT%H:%M:%S')
    my_weekday = calendar.day_name[my_datetime.weekday()] 
    my_date = my_datetime.date()
    my_time = my_datetime.time()
    
    if (my_weekday=='Friday' and my_time > datetime.time(14, 30, 0)):
        result = my_date + datetime.timedelta(days=3)
        result = result.strftime('%Y%m%d')
    elif (my_weekday=='Saturday'):
        result = 'weekend'
    elif (my_weekday=='Sunday' and my_time < datetime.time(14, 30, 0)):
        result = 'weekend'    
    elif (my_time > datetime.time(14, 30, 0)):
        result = my_date + datetime.timedelta(days=1)
        result = result.strftime('%Y%m%d')
    else:
        result = my_date
        result = result.strftime('%Y%m%d') 
    return result   

###########################################
###########################################
if __name__ == "__main__":
    opt = parse_option()
    print(opt)

    if not opt.local_topic_model:
        # This file is the concatenation of all info files from /combined_info
        df = pd.read_csv(f'{opt.concatPath}/info_concatenate.csv', sep=',')

        df['date'] = df['TimeStamp_NY'].apply(oil_date)
        df = df[df['date']!='weekend']

        df.to_csv(f'{opt.concatPath}/date_fixed_article_level_measures.csv', index=False)
        
    else:
        from glob import glob
        
        info_files = glob(opt.concatPath + '/*_info.csv')
        info_files.sort()
        
        for k in [opt.rolling_index]:
            file = info_files[k]
            date_range = file[-22:-9]
            
            df = pd.read_csv(file, sep=',')

            df['date'] = df['TimeStamp_NY'].apply(oil_date)
            df = df[df['date']!='weekend']

            df.to_csv(f'{opt.concatPath}/{date_range}_date_fixed_article_level_measures.csv', index=False)
