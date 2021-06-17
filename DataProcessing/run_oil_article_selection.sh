#!/bin/bash

FS='/user/hw2676/code/Energy/DataProcessing/oil_RTRS'


for f in {1996..2020};
do 
	for g in {01..12};
		do  sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/oil_article_selection.py $f$g";
	done;
done
 