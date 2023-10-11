#!/bin/bash

chmod 700 topic_allocation.py
# 5G is not enough
# sge_run --grid_mem=20G --grid_ncpus=1 --grid_submit=batch "./topic_allocation.py --inputWordsPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_clustering_C --outputPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_topic_allocation --local_topic_model=True --rolling_index=266"
for rolling_index in $(seq 0 266); do
  sge_run --grid_mem=20G --grid_ncpus=1 --grid_submit=batch "./topic_allocation.py --inputWordsPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_clustering_C --outputPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_topic_allocation --local_topic_model=True --rolling_index=${rolling_index}"
done
