#!/apps/anaconda2/bin/python

"""
    Program            : refer to run_sent.sh 
    Function           : This code prepares the sentiment scores
"""

from __future__ import division
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
def get_clean(sample):
    sample = sample.lower()
    tokens = nltk.word_tokenize(sample)
    tokens = mark_negation(tokens)
    return tokens

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

def get_Neg(sample):
    words = [word for word in sample if word in Neg]
    return ' '.join(words)

def get_Pos(sample):
    words = [word for word in sample if word in Pos]
    return ' '.join(words)

def get_uncertain(sample):
    sample = sample.split(' ')
    words = [word for word in sample if word in Uncertain]
    return ' '.join(words)


import csv

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

    def rows(self, cutoff=2):
        """Helper function that returns rows of term-document matrix."""
        # Get master list of words that meet or exceed the cutoff frequency
        words = [word for word in self.doc_count \
          if self.doc_count[word] >= cutoff]
        # Return header
        yield words
        # Loop over rows
        for row in self.sparse:
            # Get word counts for all words in master list. If a word does
            # not appear in this document it gets a count of 0.
            data = [row.get(word, 0) for word in words]
            yield data

    def write_csv(self, filename, cutoff=2):
        """
        Write term-document matrix to a CSV file.
        filename is the name of the output file (e.g. 'mymatrix.csv').
        cutoff is an integer that specifies only words which appear in
        'cutoff' or more documents should be written out as columns in
        the matrix.
        """
        f = csv.writer(open(filename, 'wb'))
        for row in self.rows(cutoff=cutoff):
            f.writerow(row)

if __name__ == "__main__":
    FS='/user/user1/ra2826/oil_project/article_measures/sentiment'
    os.chdir(FS)
    with open('2014.txt', 'r') as f:
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
    LITIGIOUS_index = Line.index('LITIGIOUS')
    
    neg = Line[neg_index+1:pos_index]
    pos = Line[pos_index+1:uncertain_index]
    uncertain = Line[uncertain_index+1:LITIGIOUS_index]

    Neg = [l.lower() for l in neg]
    Pos = [l.lower() for l in pos]
    Uncertain = [l.lower() for l in uncertain]
        
    stemmer=stem.snowball.EnglishStemmer()
    stop = stopwords.words('english')
    stop_final = [stemmer.stem(l) for l in stop]

    inputpath = '/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS' 
    j = sys.argv[1]
    a = j[-10:-4]
    Temp=pd.read_csv(inputpath +'/'+j,sep=',')
    Temp['body_stem'] = Temp['augbod'].apply(get_clean1)
    Temp['body_negation'] = Temp['augbod'].apply(get_clean)
    Temp['body_Neg'] = Temp['body_negation'].apply(get_Neg)
    Temp['body_Pos'] = Temp['body_negation'].apply(get_Pos)
    Temp['body_total'] = Temp['body_stem'].apply(get_total)
    
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
    outputpath= '/NOBACKUP/scratch/ra2826/oil-project/sentiment'  
    df.to_csv(outputpath+'/'+ a +'sent.csv',index=False)


