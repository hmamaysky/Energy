#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## v: selected variables (top 2 or top 7)
## l: Lasso cross validation folds

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/Predictive_Power_of_Textual_measures/Parsimonious/9.0Version'

# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 7 2; do
            for l in 10; do
                sge_run --grid_mem=64G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/Forward_Model.py $w $f $v $l";
            done
        done
	done
done