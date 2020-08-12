#!/apps/anaconda2/bin/python
 
"""
    Program            : refer to run_entropy.sh 
    Function           : Input: info, 3gram, and 4gram files
                         Output: entropy files 
                         The code from github.com/hmamaysky/TR_text/tree/master/States/03_entropy.py is used here.   
"""


import nltk
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer  
stemmer=PorterStemmer()
from nltk import stem
from nltk.util import ngrams

import sys
import csv
import string
import pandas as pd
from collections import Counter
import re
import numpy as np
from os import listdir
from os.path import isfile, join



######################################################## 
# 
# File Paths 
# 
########################################################

ng3 = '/NOBACKUP/scratch/ra2826/oil-project/3gram'
ng4 = '/NOBACKUP/scratch/ra2826/oil-project/4gram'
info = '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS'
output = '/NOBACKUP/scratch/ra2826/oil-project/entropy'



SS='/'
datas = [info + SS + f for f in listdir(info) if isfile(join(info, f))]
datas.sort()
datas2=listdir(info)
datas2.sort()
ngr3 = [ng3 + SS + f for f in listdir(ng3) if isfile(join(ng3, f))]
ngr4 = [ng4 + SS + f for f in listdir(ng4) if isfile(join(ng4, f))]
ngr3.sort()
ngr4.sort()


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


def get_clean4(sample):
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
        trigrams = ngrams(cleanup,4)
        Trigrams.extend(['.'.join(lx) for lx in trigrams])

    return Trigrams 



class NGRAM(object):
    
    def __init__(self):
        # Set tokenizer to use for tokenizing new documents
        # The term document matrix is a sparse matrix represented as a
        # list of dictionaries. Each dictionary contains the word
        # counts for a document.
        self.sparse = []
        # Keep track of the number of documents containing the word.
        self.doc_count = {}

    def add_doc(self, document):
        """Add document to the term-document matrix."""
        # Split document up into list of strings
        words = document
        # Count word frequencies in this document
        word_counts = {}
        for word in words:
            word_counts[word] = word_counts.get(word, 0) + 1
        # Add word counts as new row to sparse matrix
        self.sparse.append(word_counts)
        # Add to total document count for each word
        for word in word_counts:
            self.doc_count[word] = self.doc_count.get(word, 0) + 1



if __name__ == "__main__":
    stemmer=stem.snowball.EnglishStemmer()
    stop = stopwords.words('english')
    fnum = int(sys.argv[1])
    Temp=pd.read_csv(datas[fnum],sep=',')
    Temp['gram4'] = Temp['augbod'].apply(get_clean4)
    ngram4s = NGRAM()
    for f in Temp['gram4']:
        ngram4s.add_doc(f)

    # fnum 27 indicates the number of months to use as the trial data to use for entropy calculations
    if fnum >= 27:
        stop4_dict = {}
        for l in ngr4[fnum-27:fnum-3]:
            with open(l, 'rb') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for i,j,k in spamreader:
                    if i!='':
                        stop4_dict[j] = stop4_dict.get(j,0) + int(k)


        stop3_dict = {}
        for l in ngr3[fnum-27:fnum-3]:
            with open(l, 'rb') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for i,j,k in spamreader:
                    if i!='':
                        stop3_dict[j] = stop3_dict.get(j,0) + int(k)

        entropy = []
        for a in ngram4s.sparse:
            p = [float(l)/sum(a.values()) for l in a.values()]
            M = []
            for l in a.keys():
                b = '.'.join(l.split('.')[:3])
                temp = float((stop4_dict.get(l,0) + 1))/(stop3_dict.get(b,0) + 10)
                M.append(temp)
            m = [-p[i] * np.log(M[i]) for i in range(len(M))]
            entropy.append(sum(m))
    else:
        entropy = [np.nan] * len(Temp)

    entpd = pd.DataFrame({'entropy':entropy})
    temp2 = Temp[['Id']]
    entpd = pd.concat([entpd,temp2],axis=1)
    jj=datas2[fnum]
    aa=jj[-10:-4]
    entpd.to_csv(output+SS+aa+'_entropy.csv',index=False)






