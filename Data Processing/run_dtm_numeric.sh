#!/bin/bash

DATA1='/work/hw2676/Energy/dtm_Clustering_C'
FS='/user/hw2676/code/Energy/article_measure/dtm'


 

cd $DATA1
for f in *.csv
do
	sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/dtm_numeric.py ${f}"
done	
