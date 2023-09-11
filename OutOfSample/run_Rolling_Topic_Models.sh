#!/bin/bash

chmod 700 Rolling_Topic_Models.py
for rolling_index in $(seq 0 266); do
  sge_run --grid_mem=32G --grid_ncpus=1 --grid_submit=batch --grid_quiet "./Rolling_Topic_Models.py --rolling_index=${rolling_index}"
done
