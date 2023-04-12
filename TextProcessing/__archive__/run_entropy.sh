#!/bin/bash
FS='/user/hw2676/code/Energy/DataProcessing/article_measure/entropy'

for f in {0..300}
do
	X=${f}
	sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/entropy.py ${f}"
done	
