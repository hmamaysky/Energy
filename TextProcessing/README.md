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
4.  Calculate the entropy (first step takes 6h without parallelizing; second step takes 75min)
```
chmod 700 ngram.py

grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py


chmod 700 entropy.py

./entropy.py

./sanity_check.py --check=entropy

```
5.  Calculate the sentiments and count the total number of words in each article after cleaning (20min on 64 CPUs)
```
chmod 700 sentcode.py

pip install textmining3

./sentcode.py

./sanity_check.py --check=sentiment

./sanity_check.py --check=total

```
6.  Calculates the allocation of the topics for each article
```
chmod 700 topic_allocation.py

grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py

./sanity_check.py --check=topic

```
7. Combine all the article measures and change the time to NY to create final info files
```
chmod 700 info.py

./info.py

```
All outputs are stored under `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`.

8. Concat all info files
```
chmod 700 concat.py

./concat.py

```
9. Fix the dates on info files based on the oil price eastern closing time 

```
chmod 700 date_fixed_measures.py

./date_fixed_measures.py
```
10. Aggregate from transcripts to daily measure (we take weighted average of each measure where the weights are word counts of a transcript)

```
chmod 700 agg_daily.py

./agg_daily.py
```

*** to run the code for cosine file and clustering ***

1. Process the dtm files for the cosine code 
```
chmod 700 dtm_numeric.py

./dtm_numeric.py

```
2. Prepare the cosine file
```
chmod 700 cosine.py

./cosine.py
```
3. use Louvain algorithm for clustering

*** to process texts based on two new topic allocations ***

1. Create symbolic links to old files
```
chmod 700 create_links.sh

./create_links.sh

```
2. Repeat text-processing to get daily aggregate data
```
chmod 700 repeat_new_topic.sh

./repeat_new_topic.sh

```
