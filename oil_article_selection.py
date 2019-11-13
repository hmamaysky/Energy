#!/apps/anaconda2/bin/python

"""

Created on 2019
@author: Roya Arab Loodaricheh

Program: Refer to run_oil_article_selection.sh 
Description  : 	Input is the info-file of all articles 
				First, get the articles with energy tag
				Then, in articles with similar PNAC, get the latest one
Return    : 	dataframe of selected articles' information

"""  
import pandas as pd
import sys


energyq = pd.read_csv('energytag.csv',sep=',')
energyq = energyq.energytag.tolist()
energyq = map(lambda x: 'N2:'+x.upper(),energyq)



inputpath= '/NOBACKUP/scratch/ra2826/oil-project/info'  
outputpath= '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS'  

j = sys.argv[1]
Temp = pd.read_csv(inputpath + '/' + 'RTRS-' + j + 'proc.csv', sep=',')



data1 =[]
for i in range(len(Temp)):
    subs = Temp.loc[i,'subject']
    qcheck = []
    qcheck = [elem in subs for elem in energyq]
    l = True in qcheck
    data1.append(l)
    
Temp['energyq_check'] = data1

Temp = Temp[Temp['energyq_check']==True]

Temp = Temp.sort_values('TimeStamp').groupby('PNAC')
# keep the first article in chain (before revisions)
Temp1 = Temp.first().reset_index()  
Temp1 = Temp1[Temp1['augbod'].notnull()].reset_index()

Temp1 = Temp1[['Id', 'TimeStamp', 'headline', 'subject', 'augbod']]
    
Temp1.to_csv(outputpath + '/oil_RTRS_' + j + '.csv' , encoding = 'utf-8', index=False)

