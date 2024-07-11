#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/New_Var/2.0Version/analysis_codes'

# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 1 2; do
            sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/OLS_Model.py $w $f $v";
        done
	done
done
