#!/user/kh3191/.conda/envs/nlp/bin/python
 
"""
    Function           : This code converts the words and article ids in dtm files to numbers.
                        The reason is: this is the format that the sparse matrix in cosine code gets
"""

import pandas as pd
import os
from tqdm import tqdm

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputWordsPath', type=str, 
           default='clustering_C.csv')
           #default='2018-05-04 energy word grouping 387 words.xlsx')
    parser.add_argument('--dtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C')
    parser.add_argument('--concatInfoPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv',
           help='Only used to get IDs of selected articles, without using other information, e.g., topic allocations.')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441')
           #default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


if opt.inputWordsPath.endswith('csv'):
    words_test = pd.read_csv(opt.inputWordsPath, sep=',', names=['word','Topic','freq'])
    word_set = words_test['word'][1:].tolist()
elif opt.inputWordsPath.endswith('xlsx'):
    words_test = pd.read_excel(opt.inputWordsPath)
    word_set = words_test['Word'].tolist()
    
dic_word = {k: v for v, k in enumerate(word_set)}
def func_word(item):
    return dic_word.get(item, -1)

df_selected = pd.read_csv(opt.concatInfoPath, sep=',')
id_set = df_selected.Id.tolist()
dic_id = {k: v for v, k in enumerate(id_set)}
def func_id(item):
    return dic_id.get(item, -1)


if __name__ == "__main__":
    
    for file in tqdm(os.listdir(opt.dtmPath)):

        df = pd.read_csv(f"{opt.dtmPath}/{file}")
        df = pd.melt(df, id_vars=list(df)[:2], value_vars=list(df)[2:], var_name='words', value_name='freq')
        df = df[df['freq']!=0]

        word_column = df.words.tolist()
        id_column = df.Id.tolist()
        freq_column = df.freq.tolist()

        V = freq_column
        I = list(map(func_id,id_column))
        J = list(map(func_word,word_column))

        df_out = pd.DataFrame()
        df_out['freq'] = V
        df_out['Id'] = I
        df_out['words'] = J

        df_out.to_csv(f"{opt.outputPath}/{file}",index=False)
