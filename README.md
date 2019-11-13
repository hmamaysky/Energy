# Energy

Contains codes for Energy Project


*** to run the code for info files***

1. Generate monthly csv info files (raw info files)
```
chmod u+x 

./

```
2. Select the oil articles and generate the oil-related monthly info files
```
chmod u+x oil_article_selection.py

chmod u+x run_oil_article_selection.sh

./run_oil_article_selection.sh

```
3. Prepare the dtm files
```
chmod u+x dtm.py

chmod u+x run_dtm.sh 

./run_dtm.sh 

```
4.  Calculate the entropy.
```
chmod u+x ngram.py

chmod u+x run_ngram.sh

./run_ngram.sh


chmod u+x entropy.py

chmod u+x run_entropy.sh

./run_entropy.sh
```
5.  Calculate the sentiments.
```
chmod u+x sentcode.py

chmod u+x run_sent.sh

./run_sent.sh 

```
6.  Count the total number of words in each article after cleaning.
```
chmod u+x total.py

chmod u+x run_total.sh 

./run_total.sh  

```
7.  Calculates the allocation of the topics for each article.
```
chmod u+x topic_allocation.py

chmod u+x run_topic.sh  

./run_topic.sh 

```
8.  Combine all the article measures and change the time to NY.
```
chmod u+x info.py

chmod u+x run_info.sh

./run_info.sh

```
9. Concatenate all the files at '/NOBACKUP/scratch/ra2826/oil-project/info' for the next step
10. Fix the dates on info files based on the oil price eastern closing time 

```
chmod u+x date_fixed_measures.py

sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch ./date_fixed_measures.py
```

*** to run the code for cosine file and clustering***

1. Prepare the dtm file for the cosine code 
```
chmod u+x dtm_numeric.py

chmod u+x run_dtm_numeric.sh

./run_dtm_numeric.sh

```
2. Prepare the cosine file
```
chmod u+x cosine.py

sge_run --grid_mem=128G --grid_ncpus=2 --grid_submit=batch ./cosine.py
```
3. louvain.R is used for clustering.


