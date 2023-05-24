#!/user/kh3191/.conda/envs/nlp/bin/python

"""
    Function           : This code prepares the sentiment scores
"""

from nltk.tokenize import word_tokenize
from nltk.sentiment.util import mark_negation
from nltk.corpus import stopwords
from nltk import stem
stemmer = stem.snowball.EnglishStemmer()
stop = stopwords.words('english')
stop_final = [stemmer.stem(l) for l in stop]

import numpy as np
import pandas as pd
import re
import textmining
from tqdm import tqdm

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=False)

import os

######################################################## 
# 
# File Paths 
# 
########################################################

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--sentDicPath', type=str, default='2014.txt')
    parser.add_argument('--inputPath', type=str, 
           default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info')
    parser.add_argument('--outputPath', type=str, 
           default='../../../../shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/sentiment')
    opt = parser.parse_args()
    return opt

opt = parse_option()
print(opt)


######################################################## 
# 
# Functions 
# 
########################################################
def mark_neg(sample):
    tokens = word_tokenize(sample)
    tokens = mark_negation(tokens)
    return tokens

def get_clean(sample):
    document = re.sub('[^a-z]', ' ', sample)
    cleanup = document.strip().split()
    words = [word for word in cleanup if word not in set(stop)]
    result = ' '.join(words)
    return result

def get_total(sample):
    cleanup = sample.split(' ')
    return len(cleanup)


def main():
    
    def get_Neg(sample):
        words = [word for word in sample if word in Neg]
        return ' '.join(words)

    def get_Pos(sample):
        words = [word for word in sample if word in Pos]
        return ' '.join(words)

    # def get_uncertain(sample):
    #     sample = sample.split(' ')
    #     words = [word for word in sample if word in Uncertain]
    #     return ' '.join(words)
    
    with open(opt.sentDicPath, 'r') as f:
        content = f.readlines()
        Line = []
        for l in content:
            line = l.strip()
            line = re.sub('[^A-Za-z]', ' ', line)
            line_list = line.strip().split('\n')
            Line.extend(line_list)
    
    neg_index = Line.index('NEGATIVE')
    pos_index = Line.index('POSITIVE')
    uncertain_index = Line.index('UNCERTAINTY')
    #LITIGIOUS_index = Line.index('LITIGIOUS')
    
    neg = Line[neg_index+1:pos_index]
    pos = Line[pos_index+1:uncertain_index]
    #uncertain = Line[uncertain_index+1:LITIGIOUS_index]

    Neg = [l.lower() for l in neg]
    Pos = [l.lower() for l in pos]
    #Uncertain = [l.lower() for l in uncertain]
    
    # !!crucial for efficiency!!
    Neg = set(Neg)
    Pos = set(Pos)
        
    for file in tqdm(os.listdir(opt.inputPath)):
        YYYYMM = file[-15:-9]
        Temp = pd.read_csv(f"{opt.inputPath}/{file}",sep=',',encoding = "ISO-8859-1")
        Temp['augbod'] = Temp['augbod'].str.lower()
        Temp['body_stem'] = Temp['augbod'].parallel_apply(get_clean)
        Temp['body_negation'] = Temp['augbod'].parallel_apply(mark_neg)
        Temp['body_Neg'] = Temp['body_negation'].parallel_apply(get_Neg)
        Temp['body_Pos'] = Temp['body_negation'].parallel_apply(get_Pos)
        Temp['body_total'] = Temp['body_stem'].parallel_apply(get_total)

        ngram_Neg = textmining.TermDocumentMatrix()
        for f in Temp['body_Neg']:
            ngram_Neg.add_doc(f)

        ngram_Pos = textmining.TermDocumentMatrix()
        for f in Temp['body_Pos']:
            ngram_Pos.add_doc(f)

        Final = []
        Id = list(Temp['Id'])
        for i in range(len(Id)):
            temp_pos = []
            temp_neg = []
            pos = 0
            neg = 0
            for words,num in ngram_Pos.sparse[i].items():
                pos += num
            for words,num in ngram_Neg.sparse[i].items():
                neg += num
            Final.append((Id[i],(pos-neg)/Temp['body_total'][i]))


        df = pd.DataFrame(Final,columns = ['Id','sent'])   
        df.to_csv(f"{opt.outputPath}/{YYYYMM}_sent.csv",index=False)

        
if __name__ == "__main__":

    main()
