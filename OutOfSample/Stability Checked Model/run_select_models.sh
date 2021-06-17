#!/bin/bash

## y: lookback year window

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/Predictive_Power_of_Textual_measures/Fixed_model/6.0Version'

# run the code for different years and months
for y in {1..5}; do 
    sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/select_models.py $y";
done