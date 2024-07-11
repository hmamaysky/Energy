#!/bin/bash

## Part 1 for the fixed model: 2 base vars
## w: forecasting window in week
## f: updating frequency in week
## v: selected variables "base + text"
## l: Lasso cross validation folds
## x: baseline variable

# directory of the code
FS='/user/hw2676/code/Energy/Analysis/New_Var/2.0Version/analysis_codes'

declare -a baselist=("FutRet" "xomRet" "bpRet" "rdsaRet" "StikIdx" "DOilVol" "OilVol" "DInv" "DProd" "DSpot" "tnote_10y" "DFX" "sp500Ret" "basis" "WIPI_8wk" "trend" 'BEME' 'Mom' 'VIX' 'BasMom' 'DolBeta' 'InflaBeta' 'HedgPres' 'liquidity' 'OpenInt' "RPsdf_growing" "RPsdf_rolling" "vix_diff" "ovx_diff")
              
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