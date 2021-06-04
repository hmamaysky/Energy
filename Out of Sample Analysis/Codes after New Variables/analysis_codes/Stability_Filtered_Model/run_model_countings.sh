#!/bin/bash

## y: lookback year window

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/New_Var/1.0Version/analysis_codes/Stability_Filtered_Model'

# run the code for different years and months
for y in {1..5}; do 
    sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/model_countings.py $y";
done