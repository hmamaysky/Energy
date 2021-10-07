#!/apps/anaconda2/bin/python
 
"""
    Function           : This code generate the cosine file. This takes some time to run if dtm is big	
"""

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse
import numpy as np
from os import listdir
from os.path import isfile, join
import time
import datetime
start_time = time.time()

path1 = '/NOBACKUP/scratch/ra2826/oil-project/dtm_numeric'
pathin = [path1 + '/' + f for f in listdir(path1) if isfile(join(path1, f))]

frames = []
for j in pathin:
    x = pd.read_csv(j,  delimiter=',')
    frames.append(x)

df = pd.concat(frames)
print(len(df))



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
# print("--- %s seconds ---" % (time.time() - start_time))
print(similarities.shape)
del A


words_test = pd.read_csv('/user/user1/ra2826/oil_project/article_measures/dtm/clustering_C.csv' , sep=',')
word_set = words_test.word.tolist()
print(len(word_set))

df_cosine = pd.DataFrame(data=similarities, index=word_set, columns=word_set)

now = datetime.datetime.now().strftime("%Y%m%d")
df_cosine.to_csv('/NOBACKUP/scratch/ra2826/oil-project/cosine/cosine' + now + '.csv')





print("--- %s seconds ---" % (time.time() - start_time))

