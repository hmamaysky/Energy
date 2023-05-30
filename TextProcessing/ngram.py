#!/user/kh3191/.conda/envs/nlp/bin/python
 
"""
    Program            : refer to run_ngram.sh 
    Function           : Input: info files, Output: 3gram and 4gram of one month with their frequencies
"""

import re
import nltk
from nltk.corpus import stopwords
stop = stopwords.words('english')
from nltk.stem import PorterStemmer
stemmer = PorterStemmer()
from nltk.util import ngrams

import pandas as pd

import os
from tqdm import tqdm


import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


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


if __name__ == "__main__":

    for file in tqdm(os.listdir(opt.inputPath)):
        YYYYMM = file[-15:-9]
        Temp = pd.read_csv(f'{opt.inputPath}/{file}', delimiter=',')
        Temp['gram3'] = Temp['augbod'].apply(get_clean3)
        Temp['gram4'] = Temp['augbod'].apply(get_clean4)

        gram3_dict = {}
        for word_list in Temp['gram3']:
            for word in word_list:
                gram3_dict[word] = gram3_dict.get(word, 0) + 1

        gram4_dict = {}
        for word_list in Temp['gram4']:
            for word in word_list:
                gram4_dict[word] = gram4_dict.get(word, 0) + 1

        gram3 = pd.DataFrame(gram3_dict.items(),columns = ['word','freq'])
        gram4 = pd.DataFrame(gram4_dict.items(),columns = ['word','freq'])

        gram3.to_csv(f'{opt.outputPath}/3gram/{YYYYMM}_3gram.csv')
        gram4.to_csv(f'{opt.outputPath}/4gram/{YYYYMM}_4gram.csv')
