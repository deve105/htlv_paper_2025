#--------------------------------------------------------------------
# Packages
pacman::p_load(
  tidyverse,
  tidyplots
)

#--------------------------------------------------------------------
# Load data
# Directory path
dir_path <- "./2506_mutations/raw"

# List files in the directory
txt_files <- list.files(path = dir_path, pattern = "*.txt$", full.names = TRUE)

# First, check what columns are in each file
file_structures <- lapply(txt_files, function(file) {
  tryCatch({
    # Just read the header to see column structure
    df <- read.table(file, header = TRUE, sep = "\t", nrows = 1)
    return(list(file = basename(file), columns = names(df), count = length(names(df))))
  }, error = function(e) {
    return(list(file = basename(file), error = as.character(e)))
  })
})

# Print file structures to see differences
for(fs in file_structures) {
  cat(fs$file, ": ", fs$count, " columns\n", sep="")
}


# Read all files into a list
file_contents <- lapply(txt_files, function(file) {
  tryCatch({
    # Read the file - modify read.table parameters based on your file format
    df <- read.table(file, 
                    header = TRUE,           # If your files have headers
                    sep = "\t",              # Tab-delimited (change if needed)
                    stringsAsFactors = FALSE,
                    comment.char = "",
                    )
    # Convert flg column to character in all data frames
    if("flg" %in% names(df)) {
      df$flg <- as.character(df$flg)
    }
                        
    
    # Add filename as attribute or column if needed
    df$source_file <- basename(file)
    
    return(df)
  }, error = function(e) {
    warning(paste("Error reading file:", file, "\n", e))
    return(NULL)  # Return NULL for files with errors
  })
})
head(file_contents)

for i in file_contents

combined_data <- do.call(rbind, file_contents)

combined_data <- bind_rows(file_contents)

column_names <- lapply(file_contents, names)
head(combined_data, 1)

snp_data <- read.table(
  "/Users/denriquez/Documents/GitHub/htlv_paper_2025/2506_mutations/raw/IRID024.hg38_rmdbsnp.txt", 
  header = TRUE,           # First row contains headers
  sep = "\t",              # Tab-delimited
  comment.char = "",      # Skip rows starting with #
  stringsAsFactors = FALSE # Keep strings as character vectors
)
snp_data
