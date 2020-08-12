#!/apps/anaconda2/bin/python

"""
Program   : refer to run_raw_info.sh 
Function  : parse the raw JSON file.
Return    : data frame of article information (raw info)
"""  
import pandas as pd
import json
import unicodedata
import sys
import csv


import pandas as pd
import json
import unicodedata
import sys
import csv

def get_info_article(JSON):    
    # ID
    Id = JSON['data']['id']
    
    # TIMESTAMP 
    TimeStamp = JSON['timestamps'][0]['timestamp']
    
    # PNAC
    PNAC = JSON['data']['altId']
    
    # HEADLINE
    headline = JSON['data']['headline']
    
    # URGENCY
    urgency = JSON['data']['urgency']
    
    # LANGUAGE
    lan = JSON['data']['language']
    
    # THE LIST OF SUBJECTS
    sub_list = JSON['data']['subjects']
    
    # Body text
    body = JSON['data']['body']
    if type(body) == unicode:
        # If the body is encoded in unicode
        body = unicodedata.normalize('NFKD', body).encode('ascii','ignore')
    if type(headline) == unicode:
        # If the headline is encoded in unicode.
        headline = unicodedata.normalize('NFKD', headline).encode('ascii','ignore')
    # Augment the body text with headline.
    if body != None:
        if headline != None:
            augmentedbody = headline + body
        else:
            augmentedbody = body
    else:
        augmentedbody = headline
    # Output data structure, each row contains the information of one article
    info = (Id,TimeStamp,PNAC,headline,urgency,lan,sub_list,body,augmentedbody)
    return info 



def gen_info(raw_file):
    """
        Function  : extract information from raw file.
        Return    : a dataframe with columns (ID, TIMESTAMP, PNAC, HEADLINE, URGENCY, LAN, SUBJECT, 
                    BODY, AUGMENTED BODY),
    """
    
    with open(raw_file) as json_file:
        data = json.load(json_file)

    sizeofjson = len(data['Items'])
    Data1 = [0]*sizeofjson

    # Loop over to extract information from the raw file.
    for i in range(sizeofjson):
        Data1[i] = get_info_article(data['Items'][i])
    ######
    # A list whose i-th element is the i-th English article.
    data1 = []
    # Filter articles by certain criteria.
    for l in Data1:
 
        lang = l[5]
        urgent = l[4]
        
        if (lang != 'en')|(int(urgent)<2):
            continue
        data1.append(l)
         
         
    # Output the data structure of article information.
    Temp1 = pd.DataFrame(data1, columns = ['Id', 'TimeStamp', 'PNAC', 'headline', 'urgency', 
        'lan', 'subject', 'body', 'augbod'])

    Temp1 = Temp1[Temp1['augbod'].notnull()].reset_index()
    Temp1 = Temp1.drop(['urgency','lan','body'],axis=1) 
    return Temp1
     


outputpath= '/NOBACKUP/scratch/ra2826/oil-project/raw_info' 

path1 = '/share/share1/share_mamaysky-glasserman/data/TRNewsArchive'

raw_file = sys.argv[1]
YYYYMM = raw_file[-15:-9]
Temp = gen_info(path1 + '/' + raw_file)
Temp.to_csv(outputpath + '/' + 'RTRS-' + YYYYMM + 'proc.csv', encoding = 'utf-8', index=False)
