#!/bin/bash

## Part 2 for the fixed model: 2 text vars
## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"
## l: Lasso cross validation folds
## x: baseline variable

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/New_Var/2.0Version/analysis_codes'

declare -a baselist=('artcount' 'entropy' 'sent' 'sCo' 'fCo' 'sGom' 'fGom' 'sEnv' 'fEnv' 'sEpg' 'fEpg' 'sBbl' 'fBbl' 'sRpc' 'fRpc' 'sEp' 'fEp' 'PCAfreq' 'PCAsent' 'PCAall')
              
# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 1; do
            for l in 10; do
                for x in ${baselist[@]}; do
                    sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/Fixed_Model_PartTwo.py $w $f $v $l $x";
                done
            done
        done
	done
done