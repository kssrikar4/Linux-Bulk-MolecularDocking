#!/bin/bash

log_dir="docking_results"
output_csv="best_affinities.csv"
temp_csv="temp_affinities.csv"

echo "Filename,Best Affinity (kcal/mol)" > "$temp_csv"

for logfile in "$log_dir"/*_output.txt; do
    best_affinity=$( \
        grep -E '^ +[0-9]+ ' "$logfile" | \
        awk '{print $2}' | \
        sort -n | \
        head -n 1 \
    )

    if [ -n "$best_affinity" ]; then
        fname=$(basename "$logfile")
        echo "\"$fname\",$best_affinity" >> "$temp_csv"
    fi
done

# Sort without csvsort (skip header, sort, re-add header)
{ head -n 1 "$temp_csv"; tail -n +2 "$temp_csv" | sort -t',' -k2 -n; } > "$output_csv"
rm "$temp_csv"

echo "Results saved to $output_csv"

