The file `2018-05-04 energy word grouping 387 words.xlsx` contains the initial hand-selected list of 387 energy words from which we constructed the list of 441 words ultimately used in the paper.
* The topic model was first estimated using the DTM for these 387 words, and later words were added to the existing topics from the initial topic model run.

Below is the directory tree for the project:
```
├── /shared/share_mamaysky-glasserman/energy_drivers
│   ├── 2020-11-16
│   │   └── DataProcessing
│   │       ├── info
│   │       │   └── YYYYMM_info.csv
│   │       └── topic_allocation
│   │           └── YYYYMMtopic_alloc.csv
│   └── 2023
│       └── DataProcessing
│           ├── article_measure
│           │   ├── 3gram
│           │   │   └── YYYYMM_3gram.csv
│           │   ├── 4gram
│           │   │   └── YYYYMM_4gram.csv
│           │   ├── entropy
│           │   │   └── YYYYMM_entropy.csv
│           │   ├── sentiment
│           │   │   └── YYYYMM_sent.csv
│           │   ├── topic_allocation
│           │   │   └── YYYYMM_topic_alloc.csv
│           │   └── total
│           │       └── YYYYMM_total.csv
│           ├── concat
│           │   ├── date_fixed_article_level_measures.csv
│           │   ├── dtm_concatenate.npz
│           │   └── info_concatenate.csv
│           ├── dtm_numeric
│           │   └── YYYYMM_dtm.csv
│           ├── oil_info
│           │   └── oil_YYYYMM_info.csv
│           ├── combined_info
│           │   └── YYYYMM_info.csv
│           ├── cosine
│           │   └── cosine.csv
│           ├── dtm_Clustering_C
│           │   └── YYYYMM_dtm.csv
│           └── info
│               └── YYYYMM_info.csv
│
└── /data/ThomsonReuters_NewsArchive
    ├── YYYY
    │   └── News.RTRS.YYYYMM.0214.txt
    └── count_of_articles.txt
```
