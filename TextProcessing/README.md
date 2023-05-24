# Energy

Contains codes for textual data processing in Energy Project. Please make sure the read through the contents of each code and change the directories to proper ones before running the analysis.


*** to run the code for info files***

1. Generate monthly csv info files (raw info files)
```
chmod 700 raw_info.py

./raw_info.py

```
2. Select the oil articles and generate the oil-related monthly raw info files
```
chmod 700 oil_article_selection.py

./oil_article_selection.py

```
3. Prepare the dtm files
```
chmod 700 dtm.py

grid_run --grid_mem=50G --grid_ncpus=16 --grid_submit=batch ./dtm.py --usePandas '' 

```
4.  Calculate the entropy
```
chmod 700 ngram.py

grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py


chmod 700 entropy.py

grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./entropy.py

```
5.  Calculate the sentiments
```
chmod 700 sentcode.py

pip install textmining3

./sentcode.py

```
6.  Count the total number of words in each article after cleaning
```
chmod u+x total.py

chmod u+x run_total.sh 

./run_total.sh  

```
7.  Calculates the allocation of the topics for each article
```
chmod 700 topic_allocation.py

grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py

```
8.  Combine all the article measures and change the time to NY to create final info files
```
chmod u+x info.py

chmod u+x run_info.sh

./run_info.sh

```
9. Concatenate all the files at '/work/hw2676/Energy/DataProcessing/info' for the next step

10. Fix the dates on info files based on the oil price eastern closing time 

```
chmod u+x date_fixed_measures.py

sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch ./date_fixed_measures.py
```
11. Aggregate from transcripts to daily measure (we take weighted average of each measure where the weights are word counts of a transcript)

```
chmod u+x agg_daily.py

sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch ./agg_daily.py
```

*** to run the code for cosine file and clustering***

1. Process the dtm files for the cosine code 
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
3. louvain.R is used for clustering


