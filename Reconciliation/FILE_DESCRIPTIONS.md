# File Descriptions

| File Name                              | Description                                                           |
|----------------------------------------|-----------------------------------------------------------------------|
| **2018-05-04 energy word grouping 387 words.xlsx** | This Excel file contains the initial hand-selected list of 387 energy words from which we constructed the list of 441 words ultimately used in the paper. |
| **clouds best acc.pdf**                | This PDF file pertains to the word clouds for the current topic model (:books:✅).                |
| **clouds best mod.pdf**                | This PDF file pertains to the word clouds for the current topic model (:books:📊).                     |
| **corr_check.ipynb**                   | This Jupyter Notebook checks the correlation between the measures using current data and the initial measures stored in the database, and writes the processed series out to Stata files.                     |
| **louvain_387.ipynb**                  | This Jupyter Notebook deals with the Louvain algorithm using 387 words. |
| **louvain_441.ipynb**                  | This Jupyter Notebook deals with the Louvain algorithm using 387 words, and allocates the remaining 54 words to existing clusters. |
| **missingHeadlines.ipynb**             | This Jupyter Notebook randomly reads lines from `missingHeadlines.txt`.         |
| **missingHeadlines.py**                | This Python script identifies all missing headlines and writes them to `missingHeadlines.txt`.            |
| **missingHeadlines.txt**               | This text file contains the headlines of all missing articles.            |
| **missingTags.ipynb**                  | This Jupyter Notebook analyzes the pattern of subject tags in the missing articles.              |
| **res_weighted_acc.pt**                | This PyTorch file stores modularity, accuracy, and cosine similarity of topic models from different seeds.                    |
| **sanity_check.ipynb**                 | This Jupyter Notebook is used for sanity checks.                      |
| **sanity_check.py**                    | This Python script is used for performing sanity checks.              |
