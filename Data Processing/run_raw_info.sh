#!/bin/bash

DATA1='/share/share1/share_mamaysky-glasserman/data/TRNewsArchive'
FS='/user/user1/ra2826/oil_project/raw_info'




cd $DATA1
for f in */*.txt;
do
	sge_run --grid_mem=32G --grid_ncpus=2 --grid_submit=batch --grid_quiet "${FS}/raw_info.py ${f}"
done	 
