#!/apps/anaconda2/bin/python

"""
    Program            : refer to run_dtm.sh 
    Function           : This code prepares the monthly dtm files
"""

import xml.etree.cElementTree as cElementTree
import time
import unicodedata
import nltk
import string
from nltk.corpus import stopwords
from nltk.stem import *
from nltk.stem import PorterStemmer  
stemmer=PorterStemmer()
  
import pandas as pd
import numpy as np
from nltk import stem
from multiprocessing import Pool
import re
from nltk.util import ngrams
import os
import datetime as dt
import random
from sklearn.feature_extraction.text import TfidfVectorizer
from collections import Counter
import sys
import csv
csv.field_size_limit(100000000)
from collections import OrderedDict


from os import listdir
from os.path import isfile, join




def get_clean0(sample):

    sample = sample.upper().split('\n')
    sample = ['' if (('NEWSDESK' in xx) & (len(sample)-yy <= 5 )) else xx for yy, xx in enumerate(sample)]
    sample = ' '.join(sample)
    sample = re.sub(r'\(.*?\)','',sample)
    sample = re.sub(r'<.*?>','',sample)
    
    sample = re.sub('CBOT|CHICAGO BOARD OF TRADE|CHICAGO MERCANTILE EXCHANGE|CME','',sample)

    sample = sample.lower()
    sample = re.sub('s&p', 'snp', sample)
    sample = re.sub('s & p', 'snp', sample)
    sample = re.sub("standard[\s]&[\s]poor's", 'snp', sample)
    sample = re.sub("standard[\s]and[\s]poor's", 'snp', sample)
    sample = re.sub("snp[\s]500", 'snp500', sample)
    sample = re.sub("dow[\s]jones[\s]industrial[\s]average", 'djia', sample)
    sample = re.sub("new[\s]york[\s]stock[\s]exchange", 'nyse', sample)
    sample = re.sub("london[\s]stock[\s]exchange", 'ftse', sample)
    sample = re.sub("stock[\s]exchange[\s]of[\s]hong[\s]kong", 'sehk', sample)
    sample = re.sub("australian[\s]stock[\s]exchange", 'asx', sample)
    sample = re.sub("fannie[\s]mae", 'fnma', sample)
    sample = re.sub("freddie[\s]mac", 'fdmc', sample)
    sample = re.sub("federal[\s]reserve", 'fed', sample)
    sample = re.sub("securities[\s]and[\s]exchange[\s]commission", 'sec', sample)
    sample = re.sub("chief[\s]executive[\s]officer", 'ceo', sample)
    sample = re.sub("chief[\s]financial[\s]officer", 'cfo', sample)
    sample = re.sub("chief[\s]operating[\s]officer", 'coo', sample)
    sample = re.sub("chief[\s]investment[\s]officer", 'cio', sample)
    sample = re.sub("vice[\s]president", 'vp', sample)
    sample = re.sub("international[\s]monetary[\s]fund", 'imf', sample)
    sample = re.sub('[0-9]{9}[0-9]+', ' _bn_ ', sample)
    sample = re.sub('[0-9]{6}[0-9]+', ' _mn_ ', sample)
    sample = re.sub('_n_[\s]+billion', ' _bn_ ', sample)
    sample = re.sub('_n_[\s]+million', ' _mn_ ', sample)
    sample = ' '.join(s for s in sample.split() if not any(c.isdigit() for c in s))

    document = re.sub('[^a-z]', ' ', sample)
    cleanup = document.strip().split()
    words = [word for word in cleanup if word not in set(stop)]
    result = ' '.join(words)
    return result


def get_clean1(sample):
    lx = get_clean0(sample)
    tokens = nltk.word_tokenize(lx)
    cleanup = [stemmer.stem(token) for token in tokens]

    return cleanup 


def get_clean2(sample):
    Trigrams = []
    sample = sample.upper().split('\n')
    sample = ' '.join(sample)
    sample = re.sub(r'\(.*?\)','',sample)
    sample = re.sub(r'<.*?>','',sample)
    total = sample.replace('?', '***').replace('!', '***').replace('.', '***').replace(':', '***').replace(';', '***').split('***')

    for l in total:
        lx = get_clean0(l)
        tokens = nltk.word_tokenize(lx)
        cleanup = [stemmer.stem(token) for token in tokens]
        trigrams = ngrams(cleanup,2)
        Trigrams.extend(['.'.join(lx) for lx in trigrams])

    return Trigrams 


def get_clean3(sample):
    Trigrams = []
    sample = sample.upper().split('\n')
    sample = ' '.join(sample)
    sample = re.sub(r'\(.*?\)','',sample)
    sample = re.sub(r'<.*?>','',sample)
    total = sample.replace('?', '***').replace('!', '***').replace('.', '***').replace(':', '***').replace(';', '***').split('***')

    for l in total:
        lx = get_clean0(l)
        tokens = nltk.word_tokenize(lx)
        cleanup = [stemmer.stem(token) for token in tokens]
        trigrams = ngrams(cleanup,3)
        Trigrams.extend(['.'.join(lx) for lx in trigrams])

    return Trigrams 




stemmer=stem.snowball.EnglishStemmer()
stop = stopwords.words('english')


j = sys.argv[1]
a = j[-10:-4]
inputwords= '/user/hw2676/code/Energy/DataProcessing/article_measure/dtm'
words_test = pd.read_csv(inputwords +'/'+'clustering_C.csv',sep=',')

words_test1 = words_test.word.tolist()
f1 = words_test1

words_test1 = map(lambda x: x ,words_test1)
my_dict = dict.fromkeys(words_test1,0)


inputpath = '/work/hw2676/Energy/oil_RTRS'

os.chdir(inputpath)

with open(j, 'rb') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                next(spamreader, None)
                data1 = []
                data2 = []
                data3 = []
                data4 = []
                data5 = []
                for row in spamreader:
                    my_dict = OrderedDict.fromkeys(words_test1,0)
                    l = row[4]

                    bod_txt = get_clean1(l)
                    
                    for word in bod_txt:
                        if word in my_dict:
                            my_dict[word] += 1


                    bod_txt = get_clean2(l)
                    
                    for word in bod_txt:
                        if word in my_dict:
                            my_dict[word] += 1


                    bod_txt = get_clean3(l)
                    
                    for word in bod_txt:
                        if word in my_dict:
                            my_dict[word] += 1

                    wtx = list(my_dict.values())

                    data1.append(row[0])
                    data2.append(row[1])
                    #data3.append(row[4])
                    data5.append(wtx)




                            




# Output the data structure of article information.
df1 = pd.DataFrame({'Id':data1, 'TimeStamp':data2, 'countwords':data5})
df1[f1] = pd.DataFrame(df1.countwords.values.tolist(), index= df1.index)
df1 = df1.drop(['countwords'],axis=1)
a = j[-10:-4]
outputpath= '/work/hw2676/Energy/dtm_Clustering_C'
df1.to_csv(outputpath+'/'+ a +'dtm.csv',index=False)





