#!/bin/bash

for file in *.txt; do #change it
    # Remove "clonality_", "_rmdup", and leading digits with an underscore
    new_name=$(echo "$file" | sed -E 's/abs_//; s/_clonality.*//; s/_rmdup.*//; s/^[0-9]{1,2}_//').txt
    
    # Rename the file only if the new name is different
    if [[ "$file" != "$new_name" ]]; then
        mv "$file" "$new_name"
    fi
done