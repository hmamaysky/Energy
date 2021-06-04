#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## l: Lasso cross validation folds
## r: round of shuffling 

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/New_Var/2.0Version/analysis_codes/Stability_Filtered_Model'

# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for l in 10; do
            for r in {01..30}; do
                sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/generate_model_coefs.py $w $f $l $r";
            done
        done
	done
done