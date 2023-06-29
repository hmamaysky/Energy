#!/user/kh3191/.conda/envs/nlp/bin/python
 
"""
    Function           : This code generate the cosine file. This takes some time to run if dtm is big	
"""

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse
import numpy as np

import os
from tqdm import tqdm

# TODO: add a boolean argument indicating using 387/441 words instead of using 'xlsx'/'csv' to decide
import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputWordsPath', type=str, 
           default='clustering_C.csv')
           #default='2018-05-04 energy word grouping 387 words.xlsx')
    parser.add_argument('--dtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441')
           #default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric')
    parser.add_argument('--outputDtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat')
    parser.add_argument('--outputCosinePath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/cosine')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)

if opt.inputWordsPath.endswith('csv'):
    words_test = pd.read_csv(opt.inputWordsPath, sep=',')
    word_set = words_test.word.tolist()
elif opt.inputWordsPath.endswith('xlsx'):
    words_test = pd.read_excel(opt.inputWordsPath)
    word_set = words_test.Word.tolist()
print(f'Number of words to create cosine-similarity matrix: {len(word_set)}')

pathin = [f"{opt.dtmPath}/{f}" for f in os.listdir(opt.dtmPath)]
frames = [pd.read_csv(j, delimiter=',') for j in tqdm(pathin)]
df = pd.concat(frames)
#df.to_csv(f'{opt.outputDtmPath}/dtm_numeric_concatenate.csv')
print(f"Length of df: {len(df)}")

#if opt.inputWordsPath.endswith('xlsx'):
df = df.query('words>=0')
print(f"Length of filtered df: {len(df)}")

word_column = df.words.tolist()
id_column = df.Id.tolist()
freq_column = df.freq.tolist()

del df #release memory

print(len(word_column))
print(len(id_column))
print(len(freq_column))
print(len(list(set(word_column))))
print(len(list(set(id_column))))


V = np.array(freq_column)
I = np.array(id_column)
J = np.array(word_column)

del freq_column #release memory
del id_column #release memory
del word_column #release memory


print(V.shape)
print(I.shape)
print(J.shape)

A = sparse.csr_matrix((V,(J,I)))
print(A.shape)
del V
del I
del J

similarities = cosine_similarity(A)
print(similarities.shape)
del A


df_cosine = pd.DataFrame(data=similarities, index=word_set, columns=word_set)
df_cosine.to_csv(f"{opt.outputCosinePath}/cosine.csv")
