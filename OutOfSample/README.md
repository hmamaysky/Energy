# Energy
## Out of Sample Analysis

Contains codes for OOS analysis for the Energy Project

### Prepare the rolling text variables for OOS analysis
1. Preprocessing for rolling topics
```
./NYtime_dtm.py
./dtm_numeric.py --dtmPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_Clustering_C --outputPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_numeric_441
```

2. Obtain rolling topic models using Louvain and Genetic Algorithm, and store the membership (30min at most for each rolling model)
```
./run_Rolling_Topic_Models.sh
```

3. Generate clustering_C.csv (refer to `louvain_rolling.ipynb`) after storing monthly dtm word frequencies
```
./concat.py --concat_info='' --concat_dtm='' --dtmPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/NYtime_dtm_Clustering_C --save_monthly_freq=True
```

4. Run topic allocation for each month based on backward-looking topic models
```
./run_topic_allocation.sh
```

5. Combine info
```
./run_info.sh
```
NOTE: No need to concat all files

6. Fix dates
```
./run_date_fixed_measures.sh
```
7. Daily aggregates (~25min)
```
./agg_daily.py --concatPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info --local_topic_model=True
```

