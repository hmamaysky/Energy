#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## l: Lasso cross validation folds
## r: round of shuffling 

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/Predictive_Power_of_Textual_measures/Fixed_model/6.0Version'

# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for l in 10; do
            for r in {01..30}; do
                sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/generate_model_coefs_v1.0.py $w $f $l $r";
            done
        done
	done
done