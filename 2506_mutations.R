# Load necessary packages
pacman::p_load(
     VariantAnnotation, vcfR, dplyr, tidyverse, 
    GenomicRanges, BSgenome.Hsapiens.UCSC.hg38,
    DT, ggplot2, tidyr, stringr, purrr
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
### Extracting files that start with "IRID"
irid_files <- all_files[str_detect(basename(all_files), "^IRID")]
batch1 <- map_dfr(irid_files, readvcf_pe)
batch1$tommo
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

#######################
## batch2
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
rmd_files
batch2 <- map_dfr(rmd_files, readvcf_pe2)
DT::datatable(batch2)
### HGVD = ,
       # ESP = ,
       # FREQ_1000G = ,
       #REF_COUNT = ref_tag,
       # ALT_COUNT = alt_tag,
DT::datatable(batch1)
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
mutations_annotated
######## Batch3
all_files
rmd_files2 <- all_files[str_detect(basename(all_files), "^rmdb") &
    !str_detect(basename(all_files), "IRID|ATL") & !str_detect(basename(all_files), "22_SS|23_SS|24_SS|24_1339|SS17|SS18|SS19|SS20" )]
rmd_files2
batch3 <- map_dfr(rmd_files2, readvcf_pe)
DT::datatable(batch3)

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

# BATCH4
all_files
rmd_files3 <- all_files[str_detect(basename(all_files), "^rmdb") &
    !str_detect(basename(all_files), "IRID|ATL") & str_detect(basename(all_files), "22_SS|23_SS|24_SS|24_1339|SS17|SS18|SS19|SS20" )]
rmd_files3
batch4 <- map_dfr(rmd_files3, readvcf_pe2)
DT::datatable(batch4)

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




mutations_annotated_final <- bind_rows(
    mutations_annotated,
    mutations_annotated2,
    mutations_annotated3,
    mutations_annotated4
)


str_replace_all(mutations_annotated_final$SAMPLE_ID, "hg38.*", "") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, ".*_IRID", "IRID") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, "rmdbsnp_", "") -> mutations_annotated_final$SAMPLE_ID
str_replace_all(mutations_annotated_final$SAMPLE_ID, "2._SS", "SS") -> mutations_annotated_final$SAMPLE_ID

mutations_annotated_final$SAMPLE_ID


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


### 
pacman::p_load(maftools)


maf_files <- list.files(
    "/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv/",
    pattern = "\\.maf$", full.names = TRUE
)

# Read each MAF file as a MAF object
maf_objects <- lapply(maf_files, function(f) {
  message(sprintf("Reading MAF file: %s", f))
  read.maf(maf = f)
})
maf_objects

rbindlist(maf_objects, fill = TRUE)
?read.maf# Merge all MAF objects into one
merged_maf <- maf_objects[[1]]
if (length(maf_objects) > 1) {
  for (i in 2:length(maf_objects)) {
    merged_maf <- mergeMaf(maf1 = merged_maf, maf2 = maf_objects[[i]])
  }
}

maf_list <- lapply(maf_files, function(f) {
  dt <- fread(f, sep="\t", header=TRUE)
  # Rename SYMBOL to Hugo_Symbol if present
  if ("SYMBOL" %in% names(dt)) setnames(dt, "SYMBOL", "Hugo_Symbol")
  dt
})

all_maf_df <- rbindlist(maf_list, fill=TRUE)
fwrite(all_maf_df, "merged_all_samples.maf", sep = "\t")
all_maf <- maftools::read.maf("merged_all_samples.maf")

merge_mafs(maf_files) -> merged_maf

# Read and merge all MAF files as data.tables
all_maf_df <- rbindlist(lapply(maf_files, fread), fill = TRUE)

# Write merged MAF to disk
fwrite(all_maf_df, "merged_all_samples.maf", sep = "\t")

# Read merged MAF with maftools
all_maf <- maftools::read.maf("merged_all_samples.maf")

# Now you can use all_maf with maftools functions
all_maf_df <- rbindlist(lapply(maf_files, read.maf), fill = TRUE)

for file in maf_files {
    message(sprintf("Reading MAF file: %s", file))
    read.maf(file) -> maf_data
    if (exists("all_maf")) {
        all_maf <- mergeMaf(all_maf, maf_data)
    } else {
        all_maf <- maf_data
    }
    message(sprintf("Merged MAF file: %s", file))
    return (all_maf)
}
###
mergedmaf = read.maf("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/merged_maf.maf")
getSampleSummary(mergedmaf) |> View()

SS01 = read.maf("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv/22_SS01_S238_L008__final.maf")
SS01
getFields(SS01)



library(data.table)
all_maf_df <- rbindlist(lapply(maf_files, fread), fill = TRUE)
fwrite(all_maf_df, "merged_all_samples.maf", sep = "\t")

all_maf <- maftools::read.maf("merged_all_samples.maf")


library(data.table)

read_and_merge_mafs <- function(maf_files) {
  # Read each MAF file
  maf_list <- lapply(maf_files, function(f) {
    tryCatch({
      dt <- fread(f, sep="\t", skip="#")
      message(sprintf("Read %s with %d rows", f, nrow(dt)))
      dt
    }, error = function(e) {
      message(sprintf("Error reading %s: %s", f, e$message))
      NULL
    })
  })
  # Merge all data.tables
  merged_maf <- rbindlist(maf_list, fill=TRUE)
  
  return(merged_maf)
}

# Example usage:
merged_maf <- read_and_merge_mafs(maf_files)
merged_maf


mergedmaf = read.maf(merged_maf)# merged_maf <- read_and_merge_mafs("path/to/mafs/*.maf")
