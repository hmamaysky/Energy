#!/bin/bash

DATA1='/NOBACKUP/scratch/ra2826/oil-project/dtm_Clustering_C'
FS='/user/user1/ra2826/oil_project/article_measures/dtm'


 

cd $DATA1
for f in *.csv
do
	sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/dtm_numeric.py ${f}"
done	