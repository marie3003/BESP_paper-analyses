#!/bin/bash

OUTPUT_FILE="tree_jobs.txt"
BASE_DIR="../results/run1/simulated_data"

# Empty existing file
> "$OUTPUT_FILE"

# Loop through each simulation type and subfolder
for sim_type in "$BASE_DIR"/*/; do
    sim_type_name=$(basename "$sim_type")
    for folder in "$sim_type"/*/; do
        folder_name=$(basename "$folder")
        tree_file="${folder}${folder_name}.trees"

        if [[ -f "$tree_file" ]]; then
            n_lines=$(wc -l < "$tree_file")
            for ((i=0; i<n_lines; i++)); do
                echo "${sim_type_name}/${folder_name} ${i}" >> "$OUTPUT_FILE"
            done
        else
            echo "Warning: Tree file not found: $tree_file" >&2
        fi
    done
done

echo "Created $OUTPUT_FILE"
