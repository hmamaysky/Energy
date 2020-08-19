#!/bin/bash

## Part 1 for the fixed model: 1 text and 1 base var
## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"
## l: Lasso cross validation folds
## x: baseline variable

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/Predictive_Power_of_Textual_measures/Fixed_model/5.0Version'

declare -a baselist=("FutRet" "xomRet" "bpRet" "rdsaRet" "DOilVol" "OilVol" "DInv" "DProd" "DSpot" "tnote_10y" "DFX" "sp500Ret" "basis" "WIPIyoy" "trend" "RPsdf_growing" "RPsdf_rolling" "vix_spx" "ovx_cl1")
              
# run the code for different years and months
for w in 8; do 
	for f in 1; do
        for v in 1; do
            for l in 10; do
                for x in ${baselist[@]}; do
                    sge_run --grid_mem=32G --grid_ncpus=8 --grid_submit=batch --grid_quiet "${FS}/Fixed_Model_PartThree.py $w $f $v $l $x";
                done
            done
        done
	done
done