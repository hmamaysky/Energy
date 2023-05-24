import re
from nltk.tokenize import word_tokenize
from nltk.sentiment.util import mark_negation
from nltk.corpus import stopwords
stop = stopwords.words('english')

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
