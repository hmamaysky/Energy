# Energy Data Processing

This directory contains code for processing textual data in the Energy Project. Before running, please read through the contents of each code and ensure that you change the directories (and the virtual environment) to the appropriate ones.

## :file_folder: Processing for Info Files

### Step :one:: Generate Monthly CSV Info Files (Raw Info Files)
- **Input:**
  - `/data/ThomsonReuters_NewsArchive`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`
```bash
chmod 700 raw_info.py
./raw_info.py
```

### Step :two:: Select Oil Articles and Generate Oil-Related Monthly Raw Info Files
- **Input:**
  - `./energytag.csv`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
```bash
chmod 700 oil_article_selection.py
./oil_article_selection.py
```

### Step :three:: Prepare the Document-Term Matrix (DTM) Files
- `./clustering_C.csv` contains a list of 441 words of manually-selected energy words.
- It also contains the topic labels from a previous version of the topic model (from [KC Fed RWP 20-20]([url](https://www.kansascityfed.org/research/research-working-papers/predicting-the-oil-market/))).

- **Input:**
  - `./clustering_C.csv`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
```bash
chmod 700 dtm.py
grid_run --grid_mem=50G --grid_ncpus=16 --grid_submit=batch ./dtm.py --usePandas=
```

## :file_folder: Creating the Global and Rolling Topic Models

### Step :four:.:one:: Process the DTM Files for the Global Topic Model
- **Input:**
  - `clustering_C.csv`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441`
```bash
chmod 700 dtm_numeric.py
./dtm_numeric.py
```

### Step :four:.:two:: Compute the Cosine Similarity Matrix for the Global Topic Model

- For the rolling topic models, the Cosine Similarity Matrices are computed within the loop of the Louvain algorithm, and hence do not require a separate step.
- **Input:**
  - `clustering_C.csv`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/cosine/cosine.csv`
```bash
chmod 700 cosine.py
./cosine.py
```

### Step :four:.:three:: Use Louvain Algorithm for Global and Rolling Topic Models
- This step generates 12,345 global topic models and 10 rolling topic models for each 5-year rolling window, updated monthly.
- Next, it selects the topic models (both global and rolling) that fall in the 99th percentile for highest modularity.
- Only the selected topic models are saved; the other 9,999 models are discarded.
- In addition to storing the topic model with the 99th percentile modularity, we also save the model most similar to `clustering_C.csv`.
```bash
chmod 700 run_Rolling_Topic_Models.sh
./run_Rolling_Topic_Models.sh
```

## :file_folder: Processing Texts Based on Two New Topic Allocations
- The following steps create article-level topic frequencies and topic sentiments for the selected topic models (both global and rolling) from the prior step. The article-level measures are then aggregated to a daily level.
- All the article-level measures will be stored in separate folders:
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`: an early version of the topic model named "clustering_C"
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_mod/article_measure`: the topic model at the 99th percentile for highest modularity
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_acc/article_measure`: the topic model most similar to "clustering_C"

### Step :five:.:one:: Calculate Entropy
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`
```bash
chmod 700 ngram.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py

chmod 700 entropy.py
./entropy.py
./sanity_check.py --check=entropy
```

### Step :five:.:two:: Calculate Sentiments and Total Word Count for Each Article

- **Input:**
  - `./2014.txt`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/sentiment`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/total`
```bash
chmod 700 sentcode.py
pip install textmining3
./sentcode.py
./sanity_check.py --check=sentiment
./sanity_check.py --check=total
```

### Step :five:.:three:: Calculate Topic Allocation for Each Article
- **Input:**
  - `./clustering_C.csv`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/topic_allocation`
```bash
chmod 700 topic_allocation.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py
./sanity_check.py --check=topic
```

### Step :six:: Combine Article Measures and Change Time to NY for Final Info Files
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`
```bash
chmod 700 info.py
./info.py
```

### Step :seven:: Concatenate All Info Files and Aggregate DTM Frequencies
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_freq`
```bash
chmod 700 concat.py
./concat.py
```

### Step :eight:: Fix Dates on Info Files Based on Oil Price Eastern Closing Time
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
```bash
chmod 700 date_fixed_measures.py
./date_fixed_measures.py
```

### Step :nine: Aggregate from Article-Level to Daily Measures
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat`
- **Output:**
  - `../data/NYtime_daily_level_measures_C_2023.csv`
```bash
chmod 700 agg_daily.py
./agg_daily.py
```


### Step :keycap_ten:.:one:: Create Symbolic Links to Old Files
- To avoid modifying old code and recalculating article measures that are independent of topic models (like entropy measures and article counts), we compute these measures once and store them in the following directories:
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/3gram`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/4gram`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/entropy`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/sentiment`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/total`
- We then create symbolic links from these directories to the corresponding folders under:
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_mod/article_measure`
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_acc/article_measure`
```bash
chmod 700 create_links.sh
./create_links.sh
```

### Step :keycap_ten:.:two:: Repeat Text-Processing to Get Daily Aggregate Measures
```bash
chmod 700 repeat_new_topic.sh
./repeat_new_topic.sh
```

### Step :keycap_ten:.:three:: Generate Stata Files
- **Output:**
  - `../data/transformed_data_physical_v19.2.dta`
  - `../data/transformed_data_prices_v19.2.dta`
