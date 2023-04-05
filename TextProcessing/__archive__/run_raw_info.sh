#!/bin/bash

DATA1='/shared/share_mamaysky-glasserman/data/TRNewsArchive'
FS='/user/hw2676/code/Energy/DataProcessing/raw_info'




cd $DATA1
for f in */*.txt; do
	sge_run --grid_mem=32G --grid_ncpus=2 --grid_submit=batch --grid_quiet "${FS}/raw_info.py ${f}"
done	 
