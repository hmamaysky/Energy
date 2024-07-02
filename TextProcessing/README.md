# Energy Data Processing

This directory contains code for processing textual data in the Energy Project. Before running, please read through the contents of each code and ensure that you change the directories (and the virtual environment) to the appropriate ones.

## :file_folder: Processing for Info Files

### Step :one:: Generate Monthly CSV Info Files (Raw Info Files)
- **Input:**
  - `/data/ThomsonReuters_NewsArchive`: This directory contains the raw news archive data from Thomson Reuters.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`: This directory will store the processed monthly CSV info files.
```bash
chmod 700 raw_info.py
./raw_info.py
```

### Step :two:: Select Oil Articles and Generate Oil-Related Monthly Raw Info Files
- **Input:**
  - `./energytag.csv`: A CSV file containing tags related to energy articles.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/info`: The directory with monthly CSV info files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: This directory will contain the selected oil-related monthly raw info files.
```bash
chmod 700 oil_article_selection.py
./oil_article_selection.py
```

### Step :three:: Prepare the Document-Term Matrix (DTM) Files
- `./clustering_C.csv` contains a list of 441 manually-selected energy-related words.
- It also contains the topic labels from a [previous version on KC Fed RWP 20-20]([url](https://www.kansascityfed.org/research/research-working-papers/predicting-the-oil-market/)).

- **Input:**
  - `./clustering_C.csv`: The CSV file with energy-related words and topic labels.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: The directory with oil-related monthly raw info files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`: This directory will store the Document-Term Matrix (DTM) files.
```bash
chmod 700 dtm.py
grid_run --grid_mem=50G --grid_ncpus=16 --grid_submit=batch ./dtm.py --usePandas=
```

## :file_folder: Creating the Global and Rolling Topic Models

### Step :four:.:one:: Process the DTM Files for the Global Topic Model
- **Input:**
  - `./clustering_C.csv`: The CSV file with energy-related words and topic labels.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`: The directory with Document-Term Matrix (DTM) files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441`: This directory will store the processed numeric DTM files.
```bash
chmod 700 dtm_numeric.py
./dtm_numeric.py
```

### Step :four:.:two:: Compute the Cosine Similarity Matrix for the Global Topic Model

- For the rolling topic models, the Cosine Similarity Matrices are computed within the loop of the Louvain algorithm, and hence do not require a separate step.
- **Input:**
  - `./clustering_C.csv`: The CSV file with energy-related words and topic labels.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441`: The directory with numeric DTM files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/cosine/cosine.csv`: This file will contain the Cosine Similarity Matrix.
```bash
chmod 700 cosine.py
./cosine.py
```

### Step :four:.:three:: Use Louvain Algorithm for Global and Rolling Topic Models
- This step generates 12,345 global topic models and 10 rolling topic models for each 5-year rolling window, updated monthly.
- Next, it selects the topic models (both global and rolling) that fall in the 99th percentile for highest modularity.
- Only the selected topic models are saved; the other 12,344 models are discarded.
- In addition to storing the topic model with the 99th percentile modularity, we also save the model most similar to `clustering_C.csv`.
```bash
chmod 700 run_Rolling_Topic_Models.sh
./run_Rolling_Topic_Models.sh
```

## :file_folder: Processing Texts Based on Different Topic Models
- The following steps create article-level topic frequencies and topic sentiments for the selected topic models (both global and rolling) from the prior step. The article-level measures are then aggregated to a daily level.
- All the article-level measures will be stored in separate folders:
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`: an early version of the topic model named "clustering_C"
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_mod/article_measure`: the topic model at the 99th percentile for highest modularity
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_acc/article_measure`: the topic model most similar to "clustering_C"

### Step :five:.:one:: Calculate Entropy
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: The directory with oil-related monthly raw info files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/3gram`: This directory will store 3-grams needed for entropy.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/4gram`: This directory will store 4-grams needed for entropy.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/entropy`: This directory will store the calculated entropy measures.
```bash
chmod 700 ngram.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./ngram.py

chmod 700 entropy.py
./entropy.py
./sanity_check.py --check=entropy
```

### Step :five:.:two:: Calculate Sentiments and Total Word Count for Each Article

- **Input:**
  - `./2014.txt`: A dictionary of sentiments.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: The directory with oil-related monthly raw info files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/sentiment`: This directory will store the calculated sentiment measures.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/total`: This directory will store the total word counts.
```bash
chmod 700 sentcode.py
pip install textmining3
./sentcode.py
./sanity_check.py --check=sentiment
./sanity_check.py --check=total
```

### Step :five:.:three:: Calculate Topic Allocation for Each Article
- **Input:**
  - `./clustering_C.csv`: The CSV file with energy-related words and topic labels.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: The directory with oil-related monthly raw info files.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`: The directory with Document-Term Matrix (DTM) files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/topic_allocation`: This directory will store the topic allocation measures.
```bash
chmod 700 topic_allocation.py
grid_run --grid_mem=50G --grid_ncpus=32 --grid_submit=batch ./topic_allocation.py
./sanity_check.py --check=topic
```

### Step :six:: Combine Article Measures and Change Time to NY for Final Info Files
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/oil_info`: The directory with oil-related monthly raw info files.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure`: The directory with various article measures.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`: This directory will store the combined info files.
```bash
chmod 700 info.py
./info.py
```

### Step :seven:: Concatenate All Info Files and Aggregate DTM Frequencies
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/combined_info`: The directory with combined info files.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_Clustering_C`: The directory with Document-Term Matrix (DTM) files.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv`: This CSV file will store the concatenated info.
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_freq`: This directory will store the aggregated DTM frequencies.
```bash
chmod 700 concat.py
./concat.py
```

### Step :eight:: Fix Dates on Info Files Based on Oil Price Eastern Closing Time
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/info_concatenate.csv`: The CSV file with concatenated info.
- **Output:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/date_fixed_article_level_measures.csv`: This CSV file will store the date-fixed info.
```bash
chmod 700 date_fixed_measures.py
./date_fixed_measures.py
```

### Step :nine: Aggregate from Article-Level to Daily Measures
- **Input:**
  - `/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/concat/date_fixed_article_level_measures.csv`: The file with concatenated info.
- **Output:**
  - `../data/NYtime_daily_level_measures_C_2023.csv`: This CSV file will store the daily aggregated measures.
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
