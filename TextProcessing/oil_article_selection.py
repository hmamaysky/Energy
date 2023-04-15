#!/user/kh3191/.conda/envs/nlp/bin/python

"""

Description  : 	Input is the info-file of all articles 
				First, get the articles with energy tag
				Then, in articles with similar PNAC, get the latest one
Return    : 	dataframe of selected articles' information

"""  
import pandas as pd
import os
from tqdm import tqdm

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--tagPath', type=str, 
                        default='energytag.csv')
    parser.add_argument('--inputPath', type=str, 
                        default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info')
    parser.add_argument('--outputPath', type=str, 
                        default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


def main():
    energyq = pd.read_csv(opt.tagPath, sep=',')
    energyq = energyq.energytag.tolist()
    energyq = set(map(lambda x: 'N2:'+x.upper(), energyq))

    for file in tqdm(os.listdir(opt.inputPath)):
        Temp = pd.read_csv(f'{opt.inputPath}/{file}', sep=',')    
        Temp['energyq_check'] = [energyq.intersection(eval(i)) != set() for i in Temp['subject']]
        Temp = Temp[Temp['energyq_check']]

        Temp = Temp.sort_values('TimeStamp').groupby('PNAC')
        # keep the first article in chain (before revisions)
        Temp1 = Temp.first().reset_index()  
        Temp1 = Temp1[Temp1['augbod'].notnull()].reset_index()

        Temp1 = Temp1[['Id', 'TimeStamp', 'headline', 'subject', 'augbod']]

        Temp1.to_csv(f'{opt.outputPath}/oil_{file}', encoding = 'utf-8', index=False)


    
if __name__ == '__main__':
    main()
