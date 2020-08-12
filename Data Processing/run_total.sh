#!/bin/bash

DATA1='/NOBACKUP/scratch/ra2826/oil-project/oil_RTRS'
FS='/user/user1/ra2826/oil_project/article_measures/total'


cd $DATA1
for f in *.csv
do
	# Strip path to article file.
	X=${f}
	sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/total.py ${f}"
done	
