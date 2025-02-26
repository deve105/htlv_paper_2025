#!/bin/bash

# Output file
output_file="concatenated_output.txt"

# Clear the output file if it exists
> "$output_file"

# Loop through all files
for file in *.txt; do
    # Get the filename without the path
    filename=$(basename "$file")
    
    # Skip the first line and add the filename as a new column
    tail -n +2 "$file" | awk -v fname="$filename" '{print $0 "\t" fname}' >> "$output_file"
done

echo "Concatenation complete. Output saved to $output_file."