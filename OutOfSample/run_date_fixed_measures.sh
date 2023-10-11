#!/bin/bash

chmod 700 date_fixed_measures.py
# sge_run --grid_mem=20G --grid_ncpus=1 --grid_submit=batch "./date_fixed_measures.py --concatPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info --local_topic_model=True --rolling_index=0"
for rolling_index in $(seq 0 266); do
  sge_run --grid_mem=20G --grid_ncpus=1 --grid_submit=batch "./date_fixed_measures.py --concatPath=/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info --local_topic_model=True --rolling_index=${rolling_index}"
done
