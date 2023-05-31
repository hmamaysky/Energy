#!/user/kh3191/.conda/envs/nlp/bin/python
"""
Adapted from Roya's codes
    Function           : This code fixes the dates on info files based on the oil price eastern closing time 
"""
import pandas as pd
import calendar
import datetime

###########################################
###########################################

# This file is the concatenation of all info files from /combined_info
df = pd.read_csv('/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv', sep=',')


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



df['date'] = df['TimeStamp_NY'].apply(oil_date)
df = df[df['date']!='weekend']

df.to_csv('/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/date_fixed_article_level_measures.csv', index=False)
