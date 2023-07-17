# File Descriptions

| File Name                         | Description                                                                                       |
|-----------------------------------|---------------------------------------------------------------------------------------------------|
| **2014.txt**                      | This file is a sentiment dictionary used to generate sentiment scores.                           |
| **agg_daily.py**                  | This Python script is used for aggregating article-level measures on a daily basis.              |
| **:books::checkered_flag:<br>clustering_C.csv**         | This CSV file records the initial topic model in the 2021-version paper.                 |
| **:books::white_check_mark:<br>clustering_C_acc.csv**  | This CSV file records the current topic model in the 2023-version paper, which is the model closest to the initial one (:books::checkered_flag:).|
| **:books::bar_chart:<br>clustering_C_mod.csv**  | This CSV file records the current topic model in the 2023-version paper, which has the highest modularity in the 99th percentile.|
| **concat.py**                     | This Python script is responsible for concatenating all info files.                              |
| **cosine.py**                     | This Python script calculates the cosine similarity matrix.                                       |
| **create_links.sh**               | This shell script creates symbolic links to existing files.                                      |
| **date_fixed_measures.py**        | This Python script is used to fix dates on information files based on the oil price eastern closing time. |
| **dtm.py**                        | This Python script is involved in preparing document-term matrix (DTM) files.                   |
| **dtm_numeric.py**                | This Python script converts words in DTM to integers.                                            |
| **energytag.csv**                 | This CSV file contains energy-related tags used to select energy-related articles from raw data.|
| **entropy.py**                    | This Python script is responsible for calculating entropy.                                       |
| **info.py**                       | This Python script combines article measures and changes the time zone to New York (NY) to create final info files. |
| **ngram.py**                      | This Python script counts n-grams in the articles.                                               |
| **oil_article_selection.py**      | This Python script selects oil-related articles from raw data.                                   |
| **raw_info.py**                   | This Python script generates monthly raw info files of all articles.                             |
| **repeat_new_topic.sh**           | This shell script repeats text-processing to get daily aggregate data based on current topic models (:books::white_check_mark:;:books::bar_chart:).|
| **sentcode.py**                   | This Python script is used to calculate sentiments and the total number of words in each article after cleaning. |
| **topic_allocation.py**           | This Python script calculates the topic allocations for each article.                            |
| **topic_stability.ipynb**         | This Jupyter Notebook file visualizes the stability of the initial topic model (:books::checkered_flag:) across time. |
| **utils.py**                      | This Python script contains functions used for text preprocessing, tokenization, and creating n-grams. |
