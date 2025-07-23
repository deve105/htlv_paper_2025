#### Analysis of. SNV/mutations in Human Samples of HTLV-1 Infected Individuals
#### This analysis is only from the pre-processed data from Tokyo University

# Load necessary packages
pacman::p_load(
     VariantAnnotation, vcfR, dplyr, tidyverse, 
    GenomicRanges, BSgenome.Hsapiens.UCSC.hg38, googlesheets4,
    DT, ggplot2, tidyr, stringr, purrr, maftools
)

### Reading all table files in the specified directory
all_files = list.files("/Users/denriquez/Documents/GitHub/htlv_paper_2025/2506_mutations/raw", full.names = TRUE)

### Function to read Table files and extract relevant columns
readvcf_pe = function(file_path) {
    df <- read_tsv(file_path, 
                    col_types=cols(),
                    col_names = TRUE)
    df$sample = gsub("\\.*", "", basename(file_path))
    return(df)
}
####-----------------------------------------------
#### This has 4 batches of data with different formats
# Extracting files that start with "IRID" and reading them
### Extracting files that start with "IRID"
irid_files <- all_files[str_detect(basename(all_files), "^IRID")]
batch1 <- map_dfr(irid_files, readvcf_pe)
batch1
# Renaming and annotating columns
mutations_annotated <- batch1 %>%
    rename(
        CHROM = `#chr`,
        POS = pos,
        FLAG = flg,
        DBSNP = dbsnp,
        COSMIC = cosmic,
        TOMMO = tommo,
        HGVD = hgvd,
        ESP = esp,
        FREQ_1000G = `1000G`,
        SYMBOL = symbol,
        FUNC = func,
        SYNO = syno,
        GENE_DETAIL = gene_detail,
        SAMPLE_COL = Sample,
        REF = ref,
        ALT = alt,
        QUAL = qual,
        REF_COUNT = ref_tag,
        ALT_COUNT = alt_tag,
        SAMPLE_ID = sample
    ) %>%
    mutate(
        # Calculate basic metrics
        TOTAL_DEPTH = REF_COUNT + ALT_COUNT,
        VAF = ALT_COUNT / TOTAL_DEPTH,
        
        # Convert quality to numeric
        QUAL = as.numeric(QUAL),
        
        # Add variant type
        VARIANT_TYPE = case_when(
            nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
            nchar(REF) > nchar(ALT) ~ "DEL",
            nchar(REF) < nchar(ALT) ~ "INS",
            TRUE ~ "COMPLEX"
        ),
        
        # Parse database information
        IN_DBSNP = !is.na(DBSNP) & DBSNP != ".",
        IN_COSMIC = !is.na(COSMIC) & COSMIC != ".",
        IN_TOMMO = !is.na(TOMMO) & TOMMO != ".",
        
        # Extract TOMMO frequency
        TOMMO_FREQ = as.numeric(gsub("TOMMO:", "", TOMMO)),
        
        # Extract COSMIC ID
        COSMIC_ID = str_extract(COSMIC, "COSM\\d+"),
        FLAG=as.character(FLAG)
    ) |>
        select(
        CHROM,
        POS,
        FLAG,
        DBSNP,
        COSMIC,
        TOMMO,
        SYMBOL,
        FUNC,
        SYNO,
        GENE_DETAIL,
        SAMPLE_COL,
        REF,
        ALT,
        QUAL,
        SAMPLE_ID,
        TOTAL_DEPTH,
        ALT_COUNT,
        REF_COUNT,
        VARIANT_TYPE,
        VAF)

#####------------------------------------------------
##### BATCH2
### This batch has a different format, so we need to read it differently
readvcf_pe2 = function(file_path) {
    df <- read_tsv(file_path, 
                    col_types=cols(),
                    col_names = FALSE)
    df$sample = gsub("\\.*", "", basename(file_path))
    return(df)
}

### Extracting files that start with "IRID"
rmd_files <- all_files[str_detect(basename(all_files), "^rmdb") &
    str_detect(basename(all_files), "IRID")]

batch2 <- map_dfr(rmd_files, readvcf_pe2)

mutations_annotated2 <- batch2 %>%
    rename(
        CHROM = 1,
        POS = 2,
        FLAG = 27,
        DBSNP = 17,
        COSMIC = 16,
        TOMMO = 31,
        SYMBOL = 7,
        FUNC = 6,
        SYNO = 9,
        GENE_DETAIL = 10,
        SAMPLE_COL = 23,
        REF = 4,
        ALT = 5,
        QUAL = 19,
        AD = 30,
        SAMPLE_ID = sample
    ) %>%
    mutate(
        # Extract the AD field (second field after splitting by ":")
        ADX = sapply(AD, function(x) strsplit(x, ":")[[1]][2]),
        # Split ADX into REF and ALT counts
        REF_COUNT = as.numeric(sapply(ADX, function(x) strsplit(x, ",")[[1]][1])),
        ALT_COUNT = as.numeric(sapply(ADX, function(x) strsplit(x, ",")[[1]][2]))
    ) |>
    mutate(
        # Calculate basic metrics
        TOTAL_DEPTH = REF_COUNT + ALT_COUNT,
        VAF = ALT_COUNT / TOTAL_DEPTH,
        # Convert quality to numeric
        QUAL = as.numeric(QUAL),
        # Add variant type
        VARIANT_TYPE = case_when(
            nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
            nchar(REF) > nchar(ALT) ~ "DEL",
            nchar(REF) < nchar(ALT) ~ "INS",
            TRUE ~ "COMPLEX"
        ),
        
        # Parse database information
        IN_DBSNP = !is.na(DBSNP) & DBSNP != ".",
        IN_COSMIC = !is.na(COSMIC) & COSMIC != ".",
        IN_TOMMO = !is.na(TOMMO) & TOMMO != ".",
        
        # Extract TOMMO frequency
        TOMMO_FREQ = as.numeric(gsub("TOMMO:", "", TOMMO)),
        
        # Extract COSMIC ID
        COSMIC_ID = str_extract(COSMIC, "COSM\\d+")
    ) |>
    select(
        CHROM,
        POS,
        FLAG,
        DBSNP,
        COSMIC,
        TOMMO,
        SYMBOL,
        FUNC,
        SYNO,
        GENE_DETAIL,
        SAMPLE_COL,
        REF,
        ALT,
        QUAL,
        SAMPLE_ID,
        TOTAL_DEPTH,
        ALT_COUNT,
        REF_COUNT,
        VARIANT_TYPE,
        VAF)
####-----------------------------------------------
#### BATCH3
# Extracting files that start with "rmdb" and do not contain "IRID"
rmd_files2 <- all_files[str_detect(basename(all_files), "^rmdb") &
    !str_detect(basename(all_files), "IRID|ATL") & !str_detect(basename(all_files), "22_SS|23_SS|24_SS|24_1339|SS17|SS18|SS19|SS20" )]
batch3 <- map_dfr(rmd_files2, readvcf_pe)

mutations_annotated3 <- batch3 %>%
    rename(
        CHROM = `#chr`,
        POS = pos,
        FLAG = filt,
        DBSNP = dbsnp,
        COSMIC = cosmic,
        TOMMO = tommo,
        SYMBOL = symbol,
        FUNC = func,
        SYNO = syno,
        GENE_DETAIL = 10,
        SAMPLE_COL = Sample,
        REF = ref,
        ALT = alt,
        QUAL = qual,
        REF_COUNT= ref_tag,
        ALT_COUNT = alt_tag,
        SAMPLE_ID = sample
    ) %>%
    mutate(
        # Calculate basic metrics
        TOTAL_DEPTH = REF_COUNT + ALT_COUNT,
        VAF = ALT_COUNT / TOTAL_DEPTH,
        
        # Convert quality to numeric
        QUAL = as.numeric(QUAL),
        
        # Add variant type
        VARIANT_TYPE = case_when(
            nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
            nchar(REF) > nchar(ALT) ~ "DEL",
            nchar(REF) < nchar(ALT) ~ "INS",
            TRUE ~ "COMPLEX"
        ),
        
        # Parse database information
        IN_DBSNP = !is.na(DBSNP) & DBSNP != ".",
        IN_COSMIC = !is.na(COSMIC) & COSMIC != ".",
        IN_TOMMO = !is.na(TOMMO) & TOMMO != ".",
        
        # Extract TOMMO frequency
        TOMMO_FREQ = as.numeric(gsub("TOMMO:", "", TOMMO)),
        
        # Extract COSMIC ID
        COSMIC_ID = str_extract(COSMIC, "COSM\\d+")
    ) |>
        select(
        CHROM,
        POS,
        FLAG,
        DBSNP,
        COSMIC,
        TOMMO,
        SYMBOL,
        FUNC,
        SYNO,
        GENE_DETAIL,
        SAMPLE_COL,
        REF,
        ALT,
        QUAL,
        SAMPLE_ID,
        TOTAL_DEPTH,
        ALT_COUNT,
        REF_COUNT,
        VARIANT_TYPE,
        VAF)
####-----------------------------------------------
#### BATCH4
# Extracting files that start with "rmdb" and do not contain "IRID"

rmd_files3 <- all_files[str_detect(basename(all_files), "^rmdb") &
    !str_detect(basename(all_files), "IRID|ATL") & str_detect(basename(all_files), "22_SS|23_SS|24_SS|24_1339|SS17|SS18|SS19|SS20" )]

batch4 <- map_dfr(rmd_files3, readvcf_pe2)

mutations_annotated4 <- batch4 %>%
    rename(
        CHROM = 1,
        POS = 2,
        FLAG = 27,
        DBSNP = 17,
        COSMIC = 16,
        TOMMO = 31,
        SYMBOL = 7,
        FUNC = 6,
        SYNO = 9,
        GENE_DETAIL = 10,
        SAMPLE_COL = 13,
        REF = 4,
        ALT = 5,
        QUAL = 26,
        AD=30,
        SAMPLE_ID = sample
    ) %>%
    mutate(
        # Extract the AD field (second field after splitting by ":")
        ADX = sapply(AD, function(x) strsplit(x, ":")[[1]][2]),
        # Split ADX into REF and ALT counts
        REF_COUNT = as.numeric(sapply(ADX, function(x) strsplit(x, ",")[[1]][1])),
        ALT_COUNT = as.numeric(sapply(ADX, function(x) strsplit(x, ",")[[1]][2]))
    ) |>
    mutate(
        # Calculate basic metrics
        TOTAL_DEPTH = REF_COUNT + ALT_COUNT,
        VAF = ALT_COUNT / TOTAL_DEPTH,
        
        # Convert quality to numeric
        QUAL = as.numeric(QUAL),
        
        # Add variant type
        VARIANT_TYPE = case_when(
            nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNV",
            nchar(REF) > nchar(ALT) ~ "DEL",
            nchar(REF) < nchar(ALT) ~ "INS",
            TRUE ~ "COMPLEX"
        ),
        
        # Parse database information
        IN_DBSNP = !is.na(DBSNP) & DBSNP != ".",
        IN_COSMIC = !is.na(COSMIC) & COSMIC != ".",
        IN_TOMMO = !is.na(TOMMO) & TOMMO != ".",
        
        # Extract TOMMO frequency
        TOMMO_FREQ = as.numeric(gsub("TOMMO:", "", TOMMO)),
        
        # Extract COSMIC ID
        COSMIC_ID = str_extract(COSMIC, "COSM\\d+")
    ) |>
        select(
        CHROM,
        POS,
        FLAG,
        DBSNP,
        COSMIC,
        TOMMO,
        SYMBOL,
        FUNC,
        SYNO,
        GENE_DETAIL,
        SAMPLE_COL,
        REF,
        ALT,
        QUAL,
        SAMPLE_ID,
        TOTAL_DEPTH,
        ALT_COUNT,
        REF_COUNT,
        VARIANT_TYPE,
        VAF)

####-----------------------------------------------
#### Final Merging of all batches
# Combine all annotated mutation data frames into one

mutations_annotated_final <- bind_rows(
    mutations_annotated,
    mutations_annotated2,
    mutations_annotated3,
    mutations_annotated4
)

#### Names Modification
str_replace_all(mutations_annotated_final$SAMPLE_ID, "hg38.*", "") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, ".*_IRID", "IRID") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, "rmdbsnp_", "") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, "2._SS", "SS") -> mutations_annotated_final$SAMPLE_ID

colnames(mutations_annotated_final)

colnames(mutations_annotated_final) <- str_replace_all(
    colnames(mutations_annotated_final), 
    c("SYMBOL" = "Hugo_Symbol",
    "CHROM" = "Chromosome",
    "POS" = "Start_Position",
    "REF" = "Reference_Allele",
    "ALT" = "Tumor_Seq_Allele2",
    "VARIANT_TYPE" = "Variant_Type",
    "FUNC" = "Functional_Consequence",
    "SYNO" = "Synonymous_Consequence",
    "SAMPLE_ID" = "Tumor_Sample_Barcode",
    "SYNO"="Variant_Classification")
)

high_quality_mutations <- mutations_annotated_final %>%
    filter(
        # Minimum depth
        TOTAL_DEPTH >= 20,
        
        # Minimum alternate allele count
        ALT_COUNT >= 5,
        
        # VAF thresholds
        VAF >= 0.005,  # Minimum VAF for somatic variants
        VAF <= 0.2,  # Avoid germline-like variants
        
        # Quality score (if available)
        if(all(!is.na(QUAL))) QUAL >= 30 else TRUE,
        
        # Focus on exonic variants
        FUNC == "exonic",
        
        # Focus on SNVs (remove this line to include indels)
        VARIANT_TYPE == "SNV"
    )
high_quality_mutations

read.maf(mutations_annotated_final) -> mutations_annotated_final_maf


# Summary statistics
cat("Total variants across all samples:", nrow(mutations_annotated_final), "\n")
cat("High-quality exonic SNVs:", nrow(high_quality_mutations), "\n")
cat("Number of samples:", length(unique(mutations_annotated_final$SAMPLE_ID)), "\n")
cat("Number of genes affected:", length(unique(mutations_annotated_final$SYMBOL[!is.na(high_quality_mutations$SYMBOL)])), "\n")

# Sample summary
sample_summary <- high_quality_mutations %>%
    dplyr::filter(!SAMPLE_ID %in% c("SS14", "SS20", "SS02", "SS04", "24_1339" )) |>
    group_by(SAMPLE_ID) %>%
    summarise(
        variants = n(),
        genes = n_distinct(SYMBOL, na.rm = TRUE),
        mean_vaf = round(mean(VAF, na.rm = TRUE), 3),
        mean_depth = round(mean(TOTAL_DEPTH, na.rm = TRUE), 1),
        .groups = "drop"
    ) %>%
    arrange(desc(variants))
sample_summary
# Gene summary
gene_summary <- high_quality_mutations %>%
        dplyr::filter(!SAMPLE_ID %in% c("SS14", "SS20", "SS02", "SS04", "24_1339" )) |>

    filter(!is.na(SYMBOL)) %>%
    count(SYMBOL, sort = TRUE) %>%
    head(20)

print("Top 20 most frequently mutated genes:")
print(gene_summary)

# Functional annotation summary
func_summary <- high_quality_mutations %>%
    count(FUNC, SYNO, sort = TRUE) 


print("Functional consequences:")
print(func_summary)

# Database overlap
db_summary <- high_quality_mutations %>%
    summarise(
        total = n(),
        in_dbsnp = sum(IN_DBSNP),
        in_cosmic = sum(IN_COSMIC),
        in_tommo = sum(IN_TOMMO),
        novel = sum(!IN_DBSNP & !IN_COSMIC),
        dbsnp_pct = round(sum(IN_DBSNP)/n()*100, 1),
        cosmic_pct = round(sum(IN_COSMIC)/n()*100, 1)
    )
print("Database annotation summary:")
print(db_summary)

library(tidyplots)
library(ggpubr)
# VAF distribution
p1 <- high_quality_mutations %>%
    ggplot(aes(x = VAF)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = "VAF Distribution", 
         x = "Variant Allele Frequency", y = "Count") +
    theme_pubr(base_size = 14, border = TRUE) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5)
    )
ggsave(p1, filename = "2506_VAF_ALL.png", width = 10, height = 10, dpi = 600)

p1  
p2
p1v2 = high_quality_mutations |>
    tidyplots::tidyplot(y=VAF) |>
    tidyplots::add_mean_area()

p1v2
# Mutations per sample
p2 <- sample_summary %>%
    ggplot(aes(x = reorder(SAMPLE_ID, variants), y = variants)) +
    geom_col(fill = "coral") +
    coord_flip() +
    labs(title = "Mutations per Sample", 
         x = "Sample", y = "Number of Mutations") +
    theme_pubr(base_size = 12, border = TRUE) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5)
    )
ggsave(p2, filename = "2506_mutations_per_sample.png", width = 10, height = 10, dpi = 600)


p2
sample_summary
sample_summary |>
    tidyplots::tidyplot(x=SAMPLE_ID, y = variants) |>
    tidyplots::add_mean_bar()

# 1. Get the top 30 most frequently mutated genes
top_genes <- high_quality_mutations %>%
    filter(!is.na(SYMBOL)) %>%
    count(SYMBOL, sort = TRUE) %>%
    slice_head(n = 30) %>%
    pull(SYMBOL)
- IRID011, IRID026 and IRID038 were excluded (no HTLV-1 clones detected)
- IRID039_S123 and SS04 (are the same)  confirmed
- IRID015 and SS14 (are the same) confirmed
- IRID062 and SS02 (are the same) confirmed
- IRID089 and SS20 (are the same) confirmed
- IRID015 and 032 are siblings

# 2. Prepare the matrix: rows = genes, columns = samples, values = mutation counts
heatmap_data <- high_quality_mutations %>%
    filter(SYMBOL %in% top_genes, !is.na(SAMPLE_ID)) %>%
    filter(!SAMPLE_ID %in% c("SS14", "SS20", "SS02", "SS04", "24_1339" )) |>
    count(SYMBOL, SAMPLE_ID) %>%
    tidyr::pivot_wider(names_from = SAMPLE_ID, values_from = n, values_fill = 0) %>%
    as.data.frame()

# 3. Convert to matrix and set row names
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$SYMBOL

# Calculate the total number of mutations per gene (row sums)
gene_totals <- rowSums(heatmap_matrix)
# Calculate the percentage for each gene (row)
gene_percent <- 100 * gene_totals / sum(heatmap_matrix)

# Create a data frame for row annotation
annotation_row <- data.frame(
    "Mutation % (Cohort)" = gene_percent
)
rownames(annotation_row) <- rownames(heatmap_matrix)

# 4. Plot the heatmap
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)
png("top30_mutation_heatmap.png", width = 400, height = 300, res = 1200)

devx = pheatmap::pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    cellwidth = 10,
    cellheight = 10,
    fontsize_col = 10,
    fontsize=7,
    main = "Top 20 Mutated Genes Across Samples",
    color = colorRampPalette(c("white", "#280697"))(100),
)
ggsave(devx, filename = "top30_mutation_heatmap.png", width = 12, height = 5, dpi = 600)






####-----------------------------------------------
####-----------------------------------------------
#### Merging MAF files from Mutect2 and Funcotator
#---1-Loading Packages
pacman::p_load(maftools, tidyverse, data.table)

#---2-Function to merge MAF objects with filename as Tumor_Sample_Barcode
merge_maf_objects <- function(maf_directory = ".", output_file = NULL) {
  # Get all MAF files in the directory
  maf_files <- list.files(path = maf_directory, 
                         pattern = "\\.(maf|MAF)$", 
                         full.names = TRUE)

  if (length(maf_files) == 0) {
    stop("No MAF files found in the specified directory")
  }
  cat("Found", length(maf_files), "MAF files:\n")
  cat(paste(basename(maf_files), collapse = "\n"), "\n\n")
  # Initialize list to store MAF objects
  maf_objects <- list()
  # Process each MAF file
  for (i in seq_along(maf_files)) {
    file_path <- maf_files[i]
    # Extract filename without extension for sample barcode
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    cat("Processing:", sample_name, "\n")
    # Read MAF file using read.maf
    tryCatch({
      maf_obj <- read.maf(maf = file_path, verbose = FALSE)
      # Get the data table from MAF object
      maf_data <- maf_obj@data
      # Update Tumor_Sample_Barcode with filename
      maf_data$Tumor_Sample_Barcode <- sample_name
      # Create new MAF object with updated sample names
      updated_maf <- read.maf(maf = maf_data, verbose = FALSE)
      # Store in list
      maf_objects[[sample_name]] <- updated_maf
      cat("  - Successfully processed", nrow(maf_data), "mutations\n")
    }, error = function(e) {
      cat("  - Error reading", sample_name, ":", e$message, "\n")
    })
  }
  if (length(maf_objects) == 0) {
    stop("No valid MAF files could be processed")
  }
  
  cat("\nMerging", length(maf_objects), "MAF objects...\n")
  
  # Extract data tables from all MAF objects
  all_maf_data <- lapply(maf_objects, function(maf_obj) {
    return(maf_obj@data)
  })
  
  # Combine all data tables
  merged_data <- rbindlist(all_maf_data, use.names = TRUE, fill = TRUE)
  
  # Create final merged MAF object
  cat("Creating merged MAF object...\n")
  merged_maf <- read.maf(maf = merged_data, verbose = FALSE)
  
  # Summary statistics
  cat("\nMerge Summary:\n")
  cat("Total mutations:", nrow(merged_maf@data), "\n")
  cat("Total samples:", length(unique(merged_maf@data$Tumor_Sample_Barcode)), "\n")
  cat("Total genes:", length(unique(merged_maf@data$Hugo_Symbol)), "\n")
  
  # Sample distribution
  sample_counts <- table(merged_maf@data$Tumor_Sample_Barcode)
  cat("\nMutations per sample:\n")
  print(sample_counts)
  
  # Variant classification summary
  cat("\nVariant classification summary:\n")
  print(table(merged_maf@data$Variant_Classification))
  
  # Save merged MAF object if output file specified
  if (!is.null(output_file)) {
    write.mafSummary(maf = merged_maf, basename = tools::file_path_sans_ext(output_file))
    fwrite(merged_maf@data, output_file, sep = "\t", quote = FALSE, na = "")
    cat("\nMerged MAF file written to:", output_file, "\n")
  }
  
  return(merged_maf)
}

## Silent mutations calculation
silent_maf_objects <- function(maf_directory = ".", output_file = NULL) {
  # Get all MAF files in the directory
  maf_files <- list.files(path = maf_directory, 
                         pattern = "\\.(maf|MAF)$", 
                         full.names = TRUE)

  if (length(maf_files) == 0) {
    stop("No MAF files found in the specified directory")
  }
  cat("Found", length(maf_files), "MAF files:\n")
  cat(paste(basename(maf_files), collapse = "\n"), "\n\n")
  # Initialize list to store MAF objects
  maf_objects <- list()
  # Process each MAF file
  for (i in seq_along(maf_files)) {
    file_path <- maf_files[i]
    # Extract filename without extension for sample barcode
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    cat("Processing:", sample_name, "\n")
    # Read MAF file using read.maf
    tryCatch({
      maf_obj <- read.maf(maf = file_path, verbose = FALSE)
      # Get the data table from MAF object
      maf_data <- maf_obj@maf.silent
      # Update Tumor_Sample_Barcode with filename
      maf_data$Tumor_Sample_Barcode <- sample_name
      # Store in list
      maf_objects[[sample_name]] <- maf_data
      cat("  - Successfully processed", nrow(maf_data), " silent mutations\n")
    }, error = function(e) {
      cat("  - Error reading", sample_name, ":", e$message, "\n")
    })
  }
  if (length(maf_objects) == 0) {
    stop("No valid MAF files could be processed")
  }
  
  cat("\nMerging", length(maf_objects), "MAF objects...\n")
  
  # Extract data tables from all MAF objects
  #all_maf_data <- lapply(maf_objects, function(maf_obj) {
   # return(maf_obj@maf.silent)
 #})
  
  # Combine all data tables
  merged_data <- rbindlist(maf_objects, use.names = TRUE, fill = TRUE)
  
  # Summary statistics
  cat("\nMerge Summary:\n")
  cat("Total mutations:", nrow(merged_data), "\n")
  cat("Total samples:", length(unique(merged_data$Tumor_Sample_Barcode)), "\n")
  cat("Total genes:", length(unique(merged_data$Hugo_Symbol)), "\n")
  
  # Sample distribution
  sample_counts <- table(merged_data$Tumor_Sample_Barcode)
  cat("\nSilent Mutations per sample:\n")
  print(sample_counts)
  
  # Variant classification summary
  cat("\nVariant classification summary:\n")
  print(table(merged_data$Variant_Classification))
  
  # Save merged MAF object if output file specified
  #if (!is.null(output_file)) {
  #  write.mafSummary(maf = merged_maf, basename = tools::file_path_sans_ext(output_file))
  #  fwrite(merged_maf@data, output_file, sep = "\t", quote = FALSE, na = "")
  #  cat("\nMerged MAF file written to:", output_file, "\n")
  #}
  
  return(merged_data)
}


silentmut=silent_maf_objects(maf_directory ="/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv/")
colnames(silentmut)
ggplot(aes(x=Start_Position), data=silentmut[Hugo_Symbol=="TRAF3"&Start_Position>102876870 & Start_Position<102876875,]) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = .7) +
    labs(title = "POPAF Distribution of Silent Mutations", 
         x = "POPAF", y = "Count") +
    theme_minimal(base_size = 14)

traf3_mut <- silentmut |>
    dplyr::filter(Start_Position==102876872) |>
  dplyr::filter(Hugo_Symbol == "TRAF3") |>
  dplyr::select(Start_Position, Tumor_Sample_Barcode, Variant_Classification, Variant_Type, POPAF, tumor_f, t_alt_count, t_ref_count, n_alt_count, n_ref_count, Reference_Allele, Tumor_Seq_Allele2)
colnames(traf3_mut)

silentmut |>
    dplyr::filter(str_detect(Tumor_Sample_Barcode,"IRID039|SS04")) |>
      dplyr::select(Hugo_Symbol , Start_Position, Tumor_Sample_Barcode, Variant_Classification, Variant_Type, POPAF, tumor_f, t_alt_count, t_ref_count, n_alt_count, n_ref_count, Reference_Allele, Tumor_Seq_Allele2) |>
    dplyr::filter(Hugo_Symbol == "TRAF3") |>
    head(30)
    dplyr::select(Start_Position, Tumor_Sample_Barcode, Variant_Classification, Variant_Type, POPAF, tumor_f, t_alt_count, t_ref_count, n_alt_count, n_ref_count, Reference_Allele, Tumor_Seq_Allele2) -> traf3_mut
traf3_mut
traf3_mut |>
    dplyr::filter(Start_Position>102877500 & Start_Position<102877600) |>
    dplyr::group_by(Tumor_Sample_Barcode, Start_Position) |>
    dplyr::summarise(Frequency = n(), .groups = "drop") |>
    print(n= Inf)

DT::datatable(traf3_mut)
    ggplot(aes(x=Start_Position))+
    geom_histogram(bins = 50, fill = "steelblue", alpha = .7) 
pacman::p_load(GenomicRanges)


lollipopPlot(data=traf3_mut, 
             gene = "TRAF3", 
             refSeqID = "NM_145725", 
             showMutationRate = TRUE, 
             showDomainLabel = TRUE, 
)+
  theme_minimal(base_size = 14)
# Create a GRanges object for TRAF3 mutations
traf3_gr <- GenomicRanges::GRanges(
  seqnames = traf3_mut$Chromosome,
  ranges = IRanges(start = traf3_mut$Start_Position, end = traf3_mut$Start_Position),
  strand = "*"
)

silentmut |>
    dplyr::filter(Hugo_Symbol =="TRAF3") |>
    dplyr::select(Start_Position) |>
    ggplot(aes(y=Start_Position)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = .7)

    #dplyr::group_by(Variant_Type, Variant_Classification) |>
    dplyr::summarise(Frequency = n(), .groups = "drop", meanAF=mean(POPAF, na.rm = TRUE)) 
silent_freq_table <- silentmut %>%
    dplyr::group_by(Hugo_Symbol, Variant_Classification) %>%
    dplyr::summarise(Frequency = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Variant_Classification, values_from = Frequency, values_fill = 0)

silent_freq_table <- silentmut %>%
    dplyr::group_by(Hugo_Symbol, Variant_Classification) %>%
    dplyr::summarise(Frequency = n(), .groups = "drop") %>%
    dplyr::arrange(desc(Frequency)) |>
    head(20)
silent_freq_table
DT::datatable(silent_freq_table)

tidyplots::tidyplot(silent_freq_table, x=Hugo_Symbol, y=Frequency, color=Variant_Classification) |>
    tidyplots::add_mean_bar() 

t1 = ggplot(silent_freq_table, aes(x = reorder(Hugo_Symbol, Frequency), y = Frequency, color=Variant_Classification)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Total Variants per Gene", x = "Gene", y = "Variant Count") +
    theme_minimal(base_size = 14)

# Optionally, save the plot
ggsave("total_variants_per_gene.png", t1, width = 10, height = 8, dpi = 300)
#---3- Run the function to merge MAF files

maf_all=merge_maf_objects(maf_directory ="/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv/")
maf_all

IRID040 = read.maf("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv/IRID040_S124_L006__final.maf")
head(IRID040@maf.silent)
subset_maf = maf_all@data |>
    dplyr::filter(!str_detect(Tumor_Sample_Barcode, "SS14|SS20|SS02|SS04|24_1339")) 
    #|>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

sub_maf = subsetMaf(maf=maf_all, tsb=subset_maf$Tumor_Sample_Barcode)
sub_maf@data = sub_maf@data |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

sub_maf@data = sub_maf@data |>
    dplyr::mutate(VAF=t_alt_count / (t_alt_count + t_ref_count))

#Variant allele frequcnies (Right bar plot)
sub_maf@data |>
    dplyr::filter(!is.na(VAF)) |>
    dplyr::group_by(Hugo_Symbol) |>
    dplyr::summarise('mean VAF' = mean(VAF, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(desc('mean VAF')) -> genes_vaf
genes_vaf 
head(genes_vaf)
getSampleSummary(sub_maf) |> View()
sub_maf@variant.classification.summary = sub_maf@variant.classification.summary |>
 dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

sub_maf@variant.type.summary = sub_maf@variant.type.summary |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

sub_maf@clinical.data = sub_maf@clinical.data |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

#### Anotate Clinical data

annocl = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=0#gid=0")
annocl = annocl |>
    dplyr::mutate(ID = str_replace(ID, "_S\\d{1,3}", "")) |>
    dplyr::mutate(ID = str_replace(ID, "^PE_", ""))
annot_clin = annocl |>
    dplyr::select(ID, disease_2)

#### Merge clinical data with MAF
sub_maf@clinical.data = sub_maf@clinical.data |>
    dplyr::left_join(annot_clin, by = c("Tumor_Sample_Barcode" = "ID")) |>
    dplyr::rename(Disease= disease_2) 

#### Test it
getSampleSummary(sub_maf) |> View()
getClinicalData(sub_maf)
library(grid)
png("2506_maf_barplot.png", width = 2400, height = 1600, res = 300)
mafbarplot(sub_maf, n=30, fontSize=0.7, legendfontSize = 0.7)
dev.off()

png("2506_maf_barplot2_summary.png", width = 2400, height = 1600, res = 300)
plotmafSummary(maf = sub_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, textSize = c(1,1))
dev.off()

png("2506_sigpw.png", width = 2400, height = 1600, res = 300)
oncoplot(maf = sub_maf, pathways = "sigpw", gene_mar = 8, fontSize = 0.7, topPathways = 3, collapsePathway = FALSE)
dev.off()

png("2506_smgpb.png", width = 2400, height = 1600, res = 300)
oncoplot(maf = sub_maf, pathways = "smgbp", gene_mar = 8, fontSize = 0.7, topPathways = 3, collapsePathway = FALSE)
dev.off()
colnames(sub_maf@data)
### Include clinical data in oncoplot
png("2506_oncoplot_clinical2.png", width = 2400, height = 1600, res = 200)
oncoplot(maf = sub_maf, 
         clinicalFeatures = "Disease",
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         , fontSize = 0.8,
         leftBarVlineCol = "blue",
         leftBarData = genes_vaf,
                  leftBarLims = c(0, 0.1))
dev.off()


#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
png("2506_KMT2D.png", width = 2400, height = 1600, res = 200)
lollipopPlot(
  maf = sub_maf,
  gene = 'KMT2D',
  labPosAngle = 45,
  showMutationRate = TRUE,
  AACol="Protein_Change",
  cBioPortal="TRUE",
  domainLabelSize=0.8,
  showDomainLabel=FALSE,
  repel=TRUE,
  collapsePosLabel=TRUE,
  legendTxtSize=1
)
dev.off()

#lollipop plot for NOTCH1, which is one of the most frequent mutated gene in Leukemia.
png("2506_NOTCH.png", width = 2400, height = 1600, res = 200)
lollipopPlot(
  maf = sub_maf,
  gene = 'NOTCH1',
  labPosAngle = 45,
  showMutationRate = TRUE,
  AACol="Protein_Change",
  cBioPortal="TRUE",
  domainLabelSize=0.8,
  showDomainLabel=FALSE,
  repel=TRUE,
  collapsePosLabel=TRUE,
  legendTxtSize=1
)
dev.off()
#lollipop plot for ACAN, which is one of the most frequent mutated gene in Leukemia.
png("2506_ACAN.png", width = 2400, height = 1600, res = 200)
lollipopPlot(
  maf = sub_maf,
  gene = 'ACAN',
  labPosAngle = 45,
  showMutationRate = TRUE,
  AACol="Protein_Change",
  cBioPortal="TRUE",
  domainLabelSize=0.8,
  showDomainLabel=FALSE,
  repel=FALSE,
  collapsePosLabel=TRUE,
  legendTxtSize=1
)
dev.off()
getFields(sub_maf) 

png("2506_rainfall_plot.png", width = 2400, height = 1600, res = 300)
rainfallPlot(maf = sub_maf, detectChangePoints = TRUE, pointSize = 0.4)
dev.off()
sub_maf


### Plot Variant Allele Frequency (VAF) for the top 20 variants
png("2506_vaf.png", width = 2400, height = 1600, res = 300)
plotVaf(maf = sub_maf, vafCol = "VAF", top = 20)
dev.off()

### Plot Somatic Interactions
png("2506_Somatic.png", width = 2400, height = 1600, res = 300)
somaticInteractions(maf = sub_maf, top = 25, pvalue = c(0.05, 0.1), fontSize=.6)
dev.off()

### Detecting cancer driver genes
mut = oncodrive(maf = sub_maf, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
mut
png("2506_Oncodrivers.png", width = 2400, height = 1600, res = 300)

plotOncodrive(res = mut, fdrCutOff = 0.05, useFraction = TRUE, labelSize = .4)
dev.off()
?plotOncodrive

#### Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = sub_maf, clinicalFeature = 'Disease')
fab.ce$groupwise_comparision[p_value < 0.05]


#### Heterogeneity
pacman::p_load(mclust)
het_h = inferHeterogeneity(maf = sub_maf, vafCol ="VAF")
het_h


#Requires BSgenome object
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

tnmatrix = trinucleotideMatrix(maf = sub_maf, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

plotApobecDiff(tnm = tnmatrix, maf = sub_maf, pVal = 0.2)

pacman::p_load(NMF)
tnsignature = estimateSignatures(mat = tnmatrix, nTry = 6)

sub_maf
