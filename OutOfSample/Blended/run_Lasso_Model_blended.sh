#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"
## l: Lasso cross validation folds

# directory of the code
FS='/user/sk5193/code/Energy/Analysis/New_Var/2.0Version/analysis_codes'

# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 1 2; do
            for l in 10; do
                sge_run --grid_mem=64G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/Lasso_Model_blended.py $w $f $v $l";
            done
        done
	done
done
