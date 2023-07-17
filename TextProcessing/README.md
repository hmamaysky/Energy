# Energy Data Processing

This directory contains code for processing textual data in the Energy Project. Before running the analysis, please read through the contents of each code and ensure that you change the directories and the virtual environment to the appropriate ones.

## Processing for Info Files

***Step :one:: Generate Monthly CSV Info Files (Raw Info Files)***
```bash
chmod 700 raw_info.py
./raw_info.py
```

***Step :two:: Select Oil Articles and Generate Oil-Related Monthly Raw Info Files***
```bash
chmod 700 oil_article_selection.py
./oil_article_selection.py
```

***Step :three:: Prepare the Document-Term Matrix (DTM) Files***
```bash
chmod 700 dtm.py
grid_run --grid_mem=50G --grid_ncpus=16 --grid_submit=batch ./dtm.py --usePandas '' 
```

***Step :four:: Calculate Entropy***
```bash
chmod 700 ngram.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py

chmod 700 entropy.py
./entropy.py
./sanity_check.py --check=entropy
```

***Step :five:: Calculate Sentiments and Total Word Count for Each Article***
```bash
chmod 700 sentcode.py
pip install textmining3
./sentcode.py
./sanity_check.py --check=sentiment
./sanity_check.py --check=total
```

***Step :six:: Calculate Topic Allocation for Each Article***
```bash
chmod 700 topic_allocation.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py
./sanity_check.py --check=topic
```

***Step :seven:: Combine Article Measures and Change Time to NY for Final Info Files***
```bash
chmod 700 info.py
./info.py
```
All outputs are stored under `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`.

***Step :eight:: Concatenate All Info Files***
```bash
chmod 700 concat.py
./concat.py
```

***Step :nine:: Fix Dates on Info Files Based on Oil Price Eastern Closing Time***
```bash
chmod 700 date_fixed_measures.py
./date_fixed_measures.py
```

***Step :one::zero:: Aggregate from Article-Level to Daily Measures***
```bash
chmod 700 agg_daily.py
./agg_daily.py
```

## Processing for Cosine File and Clustering

***Step :one:: Process the DTM Files for Cosine Code***
```bash
chmod 700 dtm_numeric.py
./dtm_numeric.py
```

***Step :two:: Prepare the Cosine File***
```bash
chmod 700 cosine.py
./cosine.py
```

***Step :three:: Use Louvain Algorithm for Clustering

## Processing Texts Based on Two New Topic Allocations

***Step :one:: Create Symbolic Links to Old Files***
```bash
chmod 700 create_links.sh
./create_links.sh
```

***Step :two:: Repeat Text-Processing to Get Daily Aggregate Measures***
```bash
chmod 700 repeat_new_topic.sh
./repeat_new_topic.sh
```
