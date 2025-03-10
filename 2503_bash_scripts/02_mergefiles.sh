#!/bin/bash

# Check if a folder is provided as an argument
if [[ -z "$1" ]]; then
    echo "Usage: $0 <folder>"
    exit 1
fi

# Folder containing the .txt files
input_folder="$1"

# Output file
output_file="concatenated_output.txt"

# Skipped files log
skipped_log="skipped_files.log"

# Clear the output file and skipped files log
> "$output_file"
> "$skipped_log"

# Check if the folder exists
if [[ ! -d "$input_folder" ]]; then
    echo "Error: Folder '$input_folder' does not exist."
    exit 1
fi

# Loop through all .txt files in the specified folder
for file in "$input_folder"/*.txt; do
    # Check if the file exists (in case no .txt files are found)
    if [[ ! -f "$file" ]]; then
        echo "No .txt files found in '$input_folder'."
        break
    fi
    
    # Count the number of lines in the file
    line_count=$(wc -l < "$file")
    
    # Skip the file if it is empty or has only one line
    if [[ "$line_count" -le 1 ]]; then
        echo "Skipping '$file' (empty or only one row)" >> "$skipped_log"
        continue
    fi
    
    # Get the filename without the path and remove the .txt extension
    filename=$(basename "$file" .txt)
    
    # Skip the first line and add the filename as a new column
    tail -n +2 "$file" | awk -v fname="Peru_${filename}" '{print $0 "\t" fname}' >> "$output_file"
done

echo "Concatenation complete. Output saved to $output_file. Skipped files logged to $skipped_log."