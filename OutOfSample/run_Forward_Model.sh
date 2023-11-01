#!/bin/bash

## w: forecasting window in week
## f: updating frequency in week
## v: selected variables (top 2 or top 7)
## l: Lasso cross validation folds

# run the code for different years and months
for w in 8; do 
    for f in 1; do
        for v in 7 6 5 4 3 2 1; do
            for l in 10; do
                sge_run --grid_mem=64G --grid_ncpus=8 --grid_submit=batch "./Forward_Model.py $w $f $v $l";
            done
        done
    done
done
