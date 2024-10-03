#!/bin/bash

for file in *; do
  if [[ "$file" == *:* ]]; then
    new_file="${file//:/_}"
    mv "$file" "$new_file"
    echo "Renamed '$file' to '$new_file'"
  fi
done
