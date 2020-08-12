#!/bin/bash
FS='/user/user1/ra2826/oil_project/article_measures/entropy'

for f in {0..277}
do
	X=${f}
	sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "${FS}/entropy.py ${f}"
done	
