#!/bin/bash

markers=("acc" "mod")

for marker in "${markers[@]}"
do
  echo "Processing marker: $marker"

  directory_path="/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing_$marker"

  echo "Executing topic_allocation.py"
  chmod 700 topic_allocation.py
  ./topic_allocation.py --inputWordsPath clustering_C_$marker.csv \
                        --outputPath $directory_path/article_measure/topic_allocation

  echo "Executing info.py"
  mkdir -p $directory_path/combined_info
  chmod 700 info.py
  ./info.py --measurePath $directory_path/article_measure \
            --outputPath $directory_path/combined_info
     
  echo "Executing concat.py"
  mkdir -p $directory_path/concat
  chmod 700 concat.py
  ./concat.py --combinedInfoPath $directory_path/combined_info \
              --outputPath $directory_path/concat
      
  echo "Executing date_fixed_measures.py"
  chmod 700 date_fixed_measures.py
  ./date_fixed_measures.py --concatPath $directory_path/concat
  
  echo "Executing agg_daily.py"
  chmod 700 agg_daily.py
  ./agg_daily.py --concatPath $directory_path/concat \
                 --outputPath NYtime_daily_level_measures_C_2023_$marker.csv
            
done