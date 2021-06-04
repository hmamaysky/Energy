#!/bin/bash

## Part 1 for the fixed model: 1 text and 1 base var
## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"
## l: Lasso cross validation folds
## x: baseline variable
## y: textual variable

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/Final_codes/20201113'

## declare -a baselist=("WIPImom_8wk" "RPsdf_growing" "RPsdf_rolling" "vix_spx" "ovx_cl1")
declare -a baselist=("trend")
declare -a textlist=('artcount' 'entropy' 'sent' 'sCo' 'fCo' 'sGom' 'fGom' 'sEnv' 'fEnv' 'sEpg' 'fEpg' 'sBbl' 'fBbl' 'sRpc' 'fRpc' 'sEp' 'fEp' 'PCAfreq' 'PCAsent' 'PCAall')             
# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 1; do
            for l in 10; do
                for x in ${baselist[@]}; do
                    for y in ${textlist[@]}; do
                        echo 'base=' $x ', text=' $y
                        sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/Fixed_Model_PartOne.py $w $f $v $l $x $y";
                    done
                done
            done
        done
	done
done