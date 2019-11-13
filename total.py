#!/apps/anaconda2/bin/python

"""
    Program            : refer to run_total.sh 
    Function           : This code counts the total number of words after cleaning for each article
"""

import xml.etree.cElementTree as cElementTree
import unicodedata
import nltk
import string
from nltk.corpus import stopwords
from nltk.stem import *
import numpy as np
import pandas as pd
import re
import textmining
import sys
from nltk import stem
from nltk.sentiment.util import *
from nltk.util import ngrams
from sklearn.feature_extraction.text import TfidfVectorizer
import os

######################################################## 
# 
# File Paths 
# 
########################################################



######################################################## 
# 
# Functions 
# 
########################################################

def get_clean1(sample):
    sample = sample.lower()
    document = re.sub('[^a-z]', ' ', sample)
    cleanup = document.strip().split()
    words = [word for word in cleanup if word not in set(stop)]
    result = ' '.join(words)
    return result

def get_total(sample):
    cleanup = sample.split(' ')
    return len(cleanup)

########################################################

import csv


if __name__ == "__main__":
    
    stemmer=stem.snowball.EnglishStemmer()
    stop = stopwords.words('english')

    inputpath = '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS' 
    j = sys.argv[1]
    a = j[-10:-4]
    Temp=pd.read_csv(inputpath +'/'+j,sep=',')
    Temp['body_stem'] = Temp['augbod'].apply(get_clean1)
    Temp['body_total'] = Temp['body_stem'].apply(get_total)
        
    Final = []
    Id = list(Temp['Id'])
    for i in range(len(Id)):
        Final.append((Id[i],Temp['body_total'][i]))


    df = pd.DataFrame(Final,columns = ['Id','total'])  
    outputpath= '/NOBACKUP/scratch/ra2826/oil-project/total'  
    df.to_csv(outputpath+'/'+ a +'total.csv',index=False)


