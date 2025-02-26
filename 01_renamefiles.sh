#!/bin/bash

for file in *.txt; do
    # Remove everything after the first underscore, remove "clonality_", remove "rmdup", and remove leading digits with an underscore
    new_name=$(echo "$file" | sed -E 's/clonality_//; s/_rmdup//; s/^[0-9]_//').txt
    
    # Rename the file (skip if the new name is the same as the old name)
    if [[ "$file" != "$new_name" ]]; then
        mv "$file" "$new_name"
    fi
done