# Energy Data Processing

This directory contains code for processing textual data in the Energy Project. Before running, please read through the contents of each code and ensure that you change the directories (and the virtual environment) to the appropriate ones.

## :file_folder: Processing for Info Files

**Step :one:: Generate Monthly CSV Info Files (Raw Info Files)***
- Input: `/data/ThomsonReuters_NewsArchive`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`
```bash
chmod 700 raw_info.py
./raw_info.py
```

**Step :two:: Select Oil Articles and Generate Oil-Related Monthly Raw Info Files***
- Input: `./energytag.csv`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
```bash
chmod 700 oil_article_selection.py
./oil_article_selection.py
```

**Step :three:: Prepare the Document-Term Matrix (DTM) Files***
- Input: `./clustering_C.csv`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
```bash
chmod 700 dtm.py
grid_run --grid_mem=50G --grid_ncpus=16 --grid_submit=batch ./dtm.py --usePandas=
```

**Step :four:: Calculate Entropy***
- Input: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`
```bash
chmod 700 ngram.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py

chmod 700 entropy.py
./entropy.py
./sanity_check.py --check=entropy
```

**Step :five:: Calculate Sentiments and Total Word Count for Each Article***
- Input: `./2014.txt`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/sentiment`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/total`
```bash
chmod 700 sentcode.py
pip install textmining3
./sentcode.py
./sanity_check.py --check=sentiment
./sanity_check.py --check=total
```

**Step :six:: Calculate Topic Allocation for Each Article***
- Input: `./clustering_C.csv`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/topic_allocation`
```bash
chmod 700 topic_allocation.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py
./sanity_check.py --check=topic
```

**Step :seven:: Combine Article Measures and Change Time to NY for Final Info Files***
- Input: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`
```bash
chmod 700 info.py
./info.py
```

**Step :eight:: Concatenate All Info Files and Aggregate DTM Frequencies***
- Input: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`&`/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_freq`
```bash
chmod 700 concat.py
./concat.py
```

**Step :nine:: Fix Dates on Info Files Based on Oil Price Eastern Closing Time***
- Input: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
- Output: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
```bash
chmod 700 date_fixed_measures.py
./date_fixed_measures.py
```

**Step :one::zero:: Aggregate from Article-Level to Daily Measures***
- Input: `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
- Output: `../data/NYtime_daily_level_measures_C_2023.csv`
```bash
chmod 700 agg_daily.py
./agg_daily.py
```

## :file_folder: Processing for Cosine File and Clustering

**Step :one:: Process the DTM Files for Cosine Code***
```bash
chmod 700 dtm_numeric.py
./dtm_numeric.py
```

**Step :two:: Prepare the Cosine File***
```bash
chmod 700 cosine.py
./cosine.py
```

**Step :three:: Use Louvain Algorithm for Clustering***

## :file_folder: Processing Texts Based on Two New Topic Allocations

**Step :one:: Create Symbolic Links to Old Files***
```bash
chmod 700 create_links.sh
./create_links.sh
```

**Step :two:: Repeat Text-Processing to Get Daily Aggregate Measures***
```bash
chmod 700 repeat_new_topic.sh
./repeat_new_topic.sh
```

**Step :three:: Generate Stata Files***
