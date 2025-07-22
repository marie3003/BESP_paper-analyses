#!/bin/bash

OUTPUT_FILE="tree_jobs.txt"
BASE_DIR="../results/pop_size_simulations/independent_homochronous"

# Empty existing file
> "$OUTPUT_FILE"

# Loop through each subfolder
for folder in "$BASE_DIR"/*/; do
    folder_name=$(basename "$folder")
    tree_file="${folder}${folder_name}.trees"

    if [[ -f "$tree_file" ]]; then
        n_lines=$(wc -l < "$tree_file")
        for ((i=0; i<n_lines; i++)); do
            echo "${folder_name} ${i}" >> "$OUTPUT_FILE"
        done
    else
        echo "Warning: Tree file not found: $tree_file" >&2
    fi
done

echo "Created $OUTPUT_FILE"
