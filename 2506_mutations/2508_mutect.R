#---1-Loading packages
pacman::p_load(maftools, tidyverse, data.table, tidyplots, corr)
pak::pkg_install("maftools")
library(tidyplots)
#----------------------------------------------
#---2-Function to merge MAF objects with filename as Tumor_Sample_Barcode
merge_maf_objects <- function(maf_directory = ".", output_file = NULL) {
  maf_files <- list.files(path = maf_directory, 
                         pattern = "\\.(maf|MAF)$", 
                         full.names = TRUE)
  if (length(maf_files) == 0) {
    stop("No MAF files found in the specified directory")
  }
  maf_objects <- list()
  for (i in seq_along(maf_files)) {
    file_path <- maf_files[i]
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    tryCatch({
      maf_obj <- read.maf(maf = file_path, verbose = FALSE)
      maf_obj@variants.per.sample$Tumor_Sample_Barcode = sample_name
      maf_obj@variant.type.summary$Tumor_Sample_Barcode = sample_name
      maf_obj@data$Tumor_Sample_Barcode = sample_name
      maf_obj@clinical.data$Tumor_Sample_Barcode = sample_name
      maf_obj@maf.silent$Tumor_Sample_Barcode = sample_name
      maf_obj@variant.classification.summary$Tumor_Sample_Barcode = sample_name
      maf_objects[[sample_name]] <- maf_obj
    }, error = function(e) {
      warning(paste("Error reading", sample_name, ":", e$message))
    })
  }
  if (length(maf_objects) == 0) {
    stop("No valid MAF files could be processed")
  }
  all_maf_data <- lapply(maf_objects, function(maf_obj) maf_obj)
  merge_mafs(all_maf_data) -> merged_maf
  return(merged_maf)

}
#----------------------------------------------
#---3-Reading data
#---metadata oficial
metadata = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=1412798579#gid=1412798579", sheet="Sheet3") |>
    mutate(status1 = ifelse(Status == "vivo", "alive", "dead/lost" )) |>
    mutate(PVLlog = as.numeric(log(PVL)))

#---MAF official
dataset = merge_maf_objects("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/Project_HTLV1 Peru_2025/2507_HTLVpaper/maf_htlv") 


#----------------------------------------------
#---4- Correlation analysis betweenthe number of mutations and PVL
correlation = dataset@variant.classification.summary |>
    mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
    mutate(ID=str_remove(ID, "^[0-9]+_")) |>
    mutate(total_nsnv=total-Missense_Mutation) |>
    right_join(metadata |> dplyr::select(PVLlog, Mean_HE, ID2, disease_2), by = c("ID"="ID2")) 

labelx = paste("Spearman rho =", 
                          round(cor(correlation$total, correlation$PVLlog, 
                                    method="spearman", use="complete.obs"), 2),
                          "\np =", 
                          round(cor.test(correlation$total, correlation$PVLlog, 
                                        method="spearman")$p.value, 3))


plot1 = correlation |>
    tidyplots::tidyplot(x=PVLlog, y=total, color=disease_2) |>
    tidyplots::add_data_points_beeswarm(size=1.5, preserve="total", alpha=.8) |>
    tidyplots::add(geom_smooth(
      method = "lm", 
      se = FALSE, 
      color = "black", 
      size=.5,
      inherit.aes=FALSE, aes(x=PVLlog, y=total),
      alpha = 0.4
      )) |>
    tidyplots::add_annotation_text(
      "Spearman rho = 0.24 \np = 0.046",
      x = min(correlation$PVLlog, na.rm=TRUE) + 1.5,
      y = max(correlation$total, na.rm=TRUE) - 1,
    ) |>
    tidyplots::adjust_x_axis("Log(PVL)") |>
    tidyplots::adjust_y_axis("Total Mutations per Case") |>
    tidyplots::adjust_legend_title("Last Disease \nStatus") |>
    tidyplots::adjust_colors(new_colors=c("#377EB8", "#E41A1C", "#a6611a")) |>
    tidyplots::adjust_font(face="bold") #|>
    tidyplots::save_plot("output/2508_PVL_vs_total_mutations.png", bg="transparent")

ggsave("output/2508_PVL_vs_total_mutations.tiff",  dpi=600)
#----------------------------------------------
#---5-MAF summary
tiff("2508_summary_maf_barplot.tiff", res = 600)
plotmafSummary(maf = dataset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, top = 20, titvRaw = FALSE)
dev.off()

#----------------------------------------------
#---5.1 Subset2 and MAF summary
sample2 = "SS14|SS04|SS02|SS20"
id_sample2 = dataset@clinical.data$Tumor_Sample_Barcode[!str_detect(dataset@clinical.data$Tumor_Sample_Barcode, sample2)]
dataset2 = maftools::subsetMaf(dataset, tsb=id_sample2)

png("2508_summary_maf2subset_barplot.png", width = 1800, height = 2400, res = 300)
plotmafSummary(maf = dataset2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, top = 20, titvRaw = FALSE)
dev.off()

#----------------------------------------------
#---5.2 Silent mutations
silent_db2 = dataset2@maf.silent
unique(silent_db2$Variant_Classification)
colnames(silent_db2)
correlations_silent = silent_db2 |>
    group_by(Tumor_Sample_Barcode, Variant_Classification) |>
    summarise(Frequency = n(), .groups = "drop") |>
    arrange(desc(Frequency)) |>
    pivot_wider(names_from = Variant_Classification, values_from = Frequency, values_fill = 0) |>
    mutate(Total_SNV = rowSums(across(where(is.numeric)))) |>
    mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
    mutate(ID=str_remove(ID, "^[0-9]+_")) |>
    right_join(metadata |> dplyr::select(PVLlog, Mean_HE, ID2, disease_2), by = c("ID"="ID2")) 

correlation_results <- correlations_silent %>%
  dplyr::select(Total_SNV, Silent, "3'UTR", "Intron", "5'Flank", "RNA", "5'UTR", "IGR", PVLlog, Mean_HE) |>
  corrr::correlate(method = "spearman") |>
  corrr::shave() |>
  corrr::fashion()
correlation_results

correlation_with_pvalues <- correlations_silent %>%
  dplyr::select(Total_SNV, Silent, "3'UTR", "Intron", "5'Flank", "RNA", "5'UTR", "IGR", PVLlog, Mean_HE) %>%
  corrr::correlate(method = "spearman", use = "complete.obs") %>%
  corrr::stretch() %>%
  dplyr::filter(!is.na(r)) %>%
  dplyr::mutate(
    p_value = purrr::map2_dbl(x, y, ~ {
      if (.x == .y) return(NA_real_)
      col_x <- correlations_silent[[.x]]
      col_y <- correlations_silent[[.y]]
      suppressWarnings(cor.test(col_x, col_y, method = "spearman", exact = FALSE)$p.value)
    })
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**", 
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::arrange(p_value)

#--------------------------------------
top_silent_genes <- silent_db2 |>
  dplyr::filter(Hugo_Symbol != "Unknown") |>
  dplyr::group_by(Hugo_Symbol) |>
  dplyr::summarise(Total_Mutations = n(), .groups = "drop") |>
  dplyr::arrange(desc(Total_Mutations)) |>
  dplyr::slice_head(n = 20) |>
  dplyr::pull(Hugo_Symbol)

graph_silent = silent_db2 |>
    dplyr::filter(Hugo_Symbol %in% top_silent_genes) |>
    group_by(Hugo_Symbol, Variant_Classification) |>
    summarise(Frequency = n(), .groups = "drop") |>
    arrange(desc(Frequency)) |>
    dplyr::filter(Hugo_Symbol!="Unknown") 

#-Median Silent mutations per sample
top_silent_genes_median <- silent_db2 |>
  dplyr::filter(Hugo_Symbol != "Unknown") |>
  dplyr::group_by(Hugo_Symbol, Tumor_Sample_Barcode) |>
  dplyr::summarise(Total_Mutations = n(), .groups = "drop") |>
  dplyr::group_by(Hugo_Symbol) |>
  dplyr::summarise(
    Median_Mutations = median(Total_Mutations),
    Q1 = quantile(Total_Mutations, 0.25),
    Q3 = quantile(Total_Mutations, 0.75),
    IQR = IQR(Total_Mutations),
    Mean_Mutations = mean(Total_Mutations),
    Total_Samples = n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(Median_Mutations)) |>
  dplyr::slice_head(n = 20)


## Graph for silent genes
grap2 = graph_silent |>
tidyplots::tidyplot(x=Hugo_Symbol, y=Frequency, color=Variant_Classification) |>
  tidyplots::add_barstack_absolute() |>
  tidyplots::adjust_x_axis(rotate_labels=TRUE) |>
  tidyplots::adjust_y_axis_title("Silent Mutation Frequency") |>
  tidyplots::adjust_legend_title("Variant Classification") |>
  tidyplots::adjust_x_axis_title("Gene Symbol") |>
  tidyplots::sort_x_axis_labels(.fun=sum, .reverse=TRUE) |>
  tidyplots::adjust_font(face="bold", fontsize=6) #|>
  tidyplots::save_plot("output/2508_Silent_mutations.png", bg="transparent")
grap2

ggsave("output/2508_Silent_mutations.tiff", plot=grap2, dpi=600, width=8, height=6)

remove_constant_cols <- function(df) {
  df %>%
    dplyr::select(where(~ length(unique(na.omit(.x))) > 1))
}

### TRAF3 mutations
silent_db2 |>
  dplyr::filter(Hugo_Symbol == "TRAF3") |>
  remove_constant_cols() |>
   mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
  mutate(ID=str_remove(ID, "^[0-9]+_")) |>
  left_join(metadata |> dplyr::select(PVLlog, Mean_HE, ID2, disease_2), by = c("ID"="ID2")) |>
  tidyplots::tidyplot(x=PVLlog, y=tumor_f, color=disease_2) |>
  tidyplots::add_data_points_beeswarm()


### TRAF3 Silent Mutations
silent_db2 |>
  dplyr::filter(Hugo_Symbol == "TRAF3") |>
  remove_constant_cols() |>
  group_by(Start_Position) |>
  summarise(Frequency = n(), .groups = "drop") |>
  arrange(desc(Frequency))
"""
 1      102877552 (rs1274174947) 38/67 (56%)
 2      102876872 (rs180709957) 37/67 (55%)
"""

### DUSP22 Silent mutations
silent_db2 |>
  dplyr::filter(Hugo_Symbol == "DUSP22") |>
  remove_constant_cols() |>
  group_by(Start_Position) |>
  summarise(Frequency = n(), .groups = "drop") |>
  arrange(desc(Frequency))
"""
320201 (rs55799519) #28/67 (41.79%)
350371 (rs3800260) #25/67 (37.31%)
"""

# Identify genes that correlate well with PVLlog following project conventions

GLM_htlv = silent_db2 |>
    mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
    mutate(ID=str_remove(ID, "^[0-9]+_")) |>
    remove_constant_cols() |>
    left_join(metadata |> dplyr::select(PVLlog, Mean_HE, ID2, disease_2), by = c("ID"="ID2")) |>
    select(Hugo_Symbol, tumor_f, PVLlog, Mean_HE, disease_2) 

# Identify genes that correlate well with PVLlog following project conventions
gene_pvl_correlations <- GLM_htlv %>%
  dplyr::filter(!is.na(tumor_f), !is.na(PVLlog)) %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::filter(n() >= 5) %>%  # Minimum 5 observations per gene for reliable correlation
  dplyr::summarise(
    n_samples = n(),
    correlation = cor(tumor_f, PVLlog, method = "spearman", use = "complete.obs"),
    p_value = suppressWarnings(cor.test(tumor_f, PVLlog, method = "spearman", exact = FALSE)$p.value),
    mean_tumor_f = round(mean(tumor_f, na.rm = TRUE), 3),
    sd_tumor_f = round(sd(tumor_f, na.rm = TRUE), 3),
    mean_pvl = round(mean(PVLlog, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    correlation_strength = case_when(
      abs(correlation) >= 0.7 ~ "Strong",
      abs(correlation) >= 0.5 ~ "Moderate",
      abs(correlation) >= 0.3 ~ "Weak",
      TRUE ~ "Very Weak"
    )
  ) %>%
  dplyr::arrange(desc(abs(correlation))) |>
  filter(p_value<0.05, n_samples > 30)


gene_pvl_correlations



##########
# Guidelines for Somatic Polymorphism Detection in Pre-Cancer
# Following htlv_paper_2025 project conventions

#----------------------------------------------
# 1. Quality Control and Preprocessing
#----------------------------------------------

# Define QC thresholds for ultra-deep sequencing
qc_thresholds <- list(
  min_depth = 500,           # Minimum read depth
  min_alt_reads = 10,         # Minimum alternative reads
  min_vaf = 0.01,            # Minimum VAF (1% for ultra-deep)
  max_vaf = 0.95,            # Maximum VAF (exclude germline)
  min_qual = 30,             # Minimum base quality
  min_mapq = 20              # Minimum mapping quality
)

# Apply QC filters to mutation data
apply_qc_filters <- function(maf_data, thresholds) {
  maf_data %>%
    dplyr::filter(
      DP >= 500,
      t_alt_count >= 10,
      tumor_f >= 0.01,
      tumor_f <= 0.95,
      !is.na(tumor_f)
    ) %>%
    dplyr::mutate(
      qc_pass = TRUE,
      detection_category = case_when(
        tumor_f >= 0.40 ~ "High frequency (>40%)",
        tumor_f >= 0.20 ~ "Intermediate frequency (20-40%)",
        tumor_f >= 0.05 ~ "Low frequency (5-20%)",
        TRUE ~ "Ultra-low frequency (<5%)"
      )
    )
}
colnames(silent_db2)
apply_qc_filters(silent_db2, qc_thresholds) -> silent_db2qc
silent_db2qc
#----------------------------------------------
# 2. Somatic vs Germline Classification
#----------------------------------------------

# Classify variants based on frequency and population data
classify_variants <- function(variant_data) {
  variant_data %>%
    dplyr::mutate(
      variant_classification = case_when(
        # Likely germline polymorphisms
        tumor_f >= 0.40 & !is.na(POPAF) & POPAF > 0.01 ~ "Germline_Polymorphism",
        tumor_f >= 0.40 & !is.na(dbSNP_RS) & dbSNP_RS != "" ~ "Likely_Germline",
        
        # Somatic mutations
        tumor_f < 0.40 & (is.na(POPAF) | POPAF < 0.001) ~ "Somatic_Mutation",
        tumor_f < 0.20 ~ "Low_Frequency_Somatic",
        
        # Uncertain significance
        TRUE ~ "Uncertain_Significance"
      ),
      evidence_level = case_when(
        !is.na(dbSNP_RS) & dbSNP_VLD == "TRUE" ~ "High",
        !is.na(POPAF) & POPAF > 0.001 ~ "Moderate",
        TRUE ~ "Low"
      )
    )
}

classify_variants(silent_db2qc) -> silent_db2qc_classified
silent_db2qc_classified |>
  dplyr::filter(Hugo_Symbol=="DUSP22") |>
  group_by(variant_classification) |>
  summarise(Frequency = n(), .groups = "drop") 
#----------------------------------------------
# 3. Pre-Cancer Specific Analysis
#----------------------------------------------

# Identify clonal evolution patterns
analyze_clonal_evolution <- function(mutation_data) {
  # Calculate mutation burden per sample
  mutation_burden <- mutation_data %>%
    dplyr::group_by(Tumor_Sample_Barcode) %>%
    dplyr::summarise(
      total_mutations = n(),
      high_freq_mutations = sum(tumor_f >= 0.20, na.rm = TRUE),
      low_freq_mutations = sum(tumor_f < 0.20 & tumor_f >= 0.05, na.rm = TRUE),
      ultra_low_freq = sum(tumor_f < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge with clinical data
  mutation_burden %>%
    dplyr::mutate(
      mutation_density = total_mutations / 1000,  # Per kb if using targeted panel
      clonal_diversity = (low_freq_mutations + ultra_low_freq) / total_mutations
    )
}

analyze_clonal_evolution(silent_db2qc_classified) -> clonal_analysis
clonal_analysis
#----------------------------------------------
# 4. Statistical Analysis Framework
#----------------------------------------------

# Perform comprehensive statistical analysis
perform_statistical_analysis <- function(classified_variants) {
  # 1. Frequency distribution analysis
  frequency_analysis <- classified_variants %>%
    dplyr::group_by(variant_classification, Hugo_Symbol) %>%
    dplyr::summarise(
      variant_count = n(),
      median_vaf = median(tumor_f, na.rm = TRUE),
      q1_vaf = quantile(tumor_f, 0.25, na.rm = TRUE),
      q3_vaf = quantile(tumor_f, 0.75, na.rm = TRUE),
      samples_affected = n_distinct(Tumor_Sample_Barcode),
      .groups = "drop"
    )
  
  # 2. Gene-level recurrence analysis
  recurrence_analysis <- classified_variants %>%
    dplyr::filter(variant_classification %in% c("Somatic_Mutation", "Low_Frequency_Somatic")) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(
      total_samples = n_distinct(Tumor_Sample_Barcode),
      recurrence_rate = total_samples / length(unique(classified_variants$Tumor_Sample_Barcode)),
      mean_vaf = mean(tumor_f, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(recurrence_rate >= 0.10) %>%  # At least 10% recurrence
    dplyr::arrange(desc(recurrence_rate))
  
  return(list(
    frequency_analysis = frequency_analysis,
    recurrence_analysis = recurrence_analysis
  ))
}

perform_statistical_analysis(silent_db2qc_classified) -> analysis_results
analysis_results
#----------------------------------------------
# 5. Visualization Framework
#----------------------------------------------

# Create comprehensive visualization plots
create_somatic_plots <- function(classified_variants) {
  # VAF distribution plot
  vaf_plot <- classified_variants %>%
    tidyplots::tidyplot(x = tumor_f, color = variant_classification) %>%
    tidyplots::add_histogram(bins = 50, alpha = 0.7) %>%
    #tidyplots::add_line_v(xintercept = c(0.05, 0.20, 0.40), 
                         #color = "red", linetype = "dashed") %>%
    tidyplots::adjust_x_axis("Variant Allele Frequency") %>%
    tidyplots::adjust_y_axis("Count") %>%
    tidyplots::adjust_legend_title("Variant\nClassification") %>%
    tidyplots::adjust_title("VAF Distribution with Classification Thresholds") %>%
    tidyplots::save_plot("output/somatic_vaf_distribution.png", 
                         width = 12, height = 8, dpi = 300)
  
  # Mutation burden by sample
  burden_plot <- classified_variants %>%
    dplyr::group_by(Tumor_Sample_Barcode, variant_classification) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    tidyplots::tidyplot(x = Tumor_Sample_Barcode, y = count, 
                        color = variant_classification) %>%
    tidyplots::add_barstack_absolute() %>%
    tidyplots::adjust_x_axis("Sample", rotate_labels = TRUE) %>%
    tidyplots::adjust_y_axis("Mutation Count") %>%
    tidyplots::adjust_legend_title("Variant\nClassification") %>%
    tidyplots::save_plot("output/mutation_burden_by_sample.png", 
                         width = 14, height = 8, dpi = 300)
  
  return(list(vaf_plot = vaf_plot, burden_plot = burden_plot))
}
create_somatic_plots(silent_db2qc_classified)
#----------------------------------------------
# 6. Reporting Framework
#----------------------------------------------

# Generate comprehensive report
generate_somatic_report <- function(classified_variants, analysis_results) {
  # Summary statistics
  summary_stats <- classified_variants %>%
    dplyr::group_by(variant_classification) %>%
    dplyr::summarise(
      total_variants = n(),
      unique_genes = n_distinct(Hugo_Symbol),
      samples_affected = n_distinct(Tumor_Sample_Barcode),
      median_vaf = round(median(tumor_f, na.rm = TRUE), 3),
      .groups = "drop"
    )
  
  # High-confidence somatic mutations
  high_confidence_somatic <- classified_variants %>%
    dplyr::filter(
      variant_classification == "Somatic_Mutation",
      tumor_f >= 0.05,
      t_alt_count >= 10
    ) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(
      mutation_count = n(),
      samples_affected = n_distinct(Tumor_Sample_Barcode),
      mean_vaf = round(mean(tumor_f), 3),
      median_vaf = round(median(tumor_f), 3),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(samples_affected), desc(mean_vaf))
  
  # Export results
  readr::write_csv(summary_stats, "output/somatic_polymorphism_summary.csv")
  readr::write_csv(high_confidence_somatic, "output/high_confidence_somatic_mutations.csv")
  readr::write_csv(classified_variants, "output/classified_variants_complete.csv")
  
  return(list(
    summary_stats = summary_stats,
    high_confidence_somatic = high_confidence_somatic
  ))
}
generate_somatic_report(silent_db2qc_classified, analysis_results) -> somatic_report
somatic_report
#----------------------------------------------
# 7. Application to HTLV Data
#----------------------------------------------

# Apply framework to current HTLV dataset
htlv_somatic_analysis <- function() {
  # Apply QC filters
  qc_filtered_data <- apply_qc_filters(silent_db2, qc_thresholds)
  
  # Classify variants
  classified_htlv <- classify_variants(qc_filtered_data)
  
  # Perform statistical analysis
  htlv_stats <- perform_statistical_analysis(classified_htlv, metadata)
  
  # Create visualizations
  htlv_plots <- create_somatic_plots(classified_htlv, metadata)
  
  # Generate report
  htlv_report <- generate_somatic_report(classified_htlv, htlv_stats, metadata)
  
  # Display key findings
  cat("=== HTLV SOMATIC POLYMORPHISM ANALYSIS ===\n")
  print(htlv_report$summary_stats)
  
  cat("\n=== HIGH-CONFIDENCE SOMATIC MUTATIONS ===\n")
  print(htlv_report$high_confidence_somatic %>% head(10))
  
  return(list(
    classified_data = classified_htlv,
    statistics = htlv_stats,
    plots = htlv_plots,
    report = htlv_report
  ))
}

#----------------------------------------------
# 8. Quality Metrics and Validation
#----------------------------------------------

# Calculate quality metrics for validation
calculate_quality_metrics <- function(variant_data) {
  quality_metrics <- variant_data %>%
    dplyr::summarise(
      total_variants = n(),
      mean_depth = mean(t_depth, na.rm = TRUE),
      mean_alt_reads = mean(t_alt_count, na.rm = TRUE),
      variants_with_dbsnp = sum(!is.na(dbSNP_RS) & dbSNP_RS != "", na.rm = TRUE),
      dbsnp_rate = variants_with_dbsnp / total_variants,
      .groups = "drop"
    )
  
  readr::write_csv(quality_metrics, "output/quality_metrics.csv")
  return(quality_metrics)
}
colnames(silent_db2)
silentgenes_csv = silent_db2 |>
  filter(Hugo_Symbol!="Unknown") |>
  group_by(Hugo_Symbol, Chromosome, Start_Position, Variant_Type, POPAF, dbSNP_ID, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2) |>
  summarise(Frequency = n()/67*100, .groups = "drop",
            mean_VAF = mean(tumor_f, na.rm = TRUE)) |>
  filter(Frequency>10) |>
  arrange(desc(Frequency))
silentgenes_csv


# Convert silentgenes_csv to publication-ready table following project conventions
pacman::p_load(flextable)
pacman::p_load(officer)

# Prepare the data for publication
publication_table <- silentgenes_csv %>%
  dplyr::mutate(
    # Format frequency as percentage with 1 decimal place
    Frequency_formatted = paste0(round(Frequency, 1), "%"),
    # Format mean VAF with 3 decimal places
    mean_VAF_formatted = round(mean_VAF, 3),
    # Create genomic position string
    Genomic_Position = paste0(Chromosome, ":", Start_Position),
    # Format allele change
    Allele_Change = paste0(Reference_Allele, ">", Tumor_Seq_Allele2),
    # Format POPAF
    POPAF_formatted = ifelse(is.na(POPAF), "Not available", 
                            ifelse(POPAF < 0.001, "<0.001", round(POPAF, 3)))
  ) %>%
  dplyr::select(
    Hugo_Symbol,
    Genomic_Position,
    Variant_Type,
    Allele_Change,
    Frequency_formatted,
    mean_VAF_formatted,
    POPAF_formatted,
    dbSNP_ID
  ) %>%
  dplyr::rename(
    "Gene Symbol" = Hugo_Symbol,
    "Genomic Position" = Genomic_Position,
    "Variant Type" = Variant_Type,
    "Allele Change" = Allele_Change,
    "Frequency (%)" = Frequency_formatted,
    "Mean VAF" = mean_VAF_formatted,
    "Population AF" = POPAF_formatted,
    "dbSNP ID" = dbSNP_ID
  )

# Create formatted flextable for publication
silent_genes_flextable <- publication_table %>%
  flextable() %>%
  # Set table style
  theme_vanilla() %>%
  # Format header
  bold(part = "header") %>%
  bg(part = "header", bg = "#E8E8E8") %>%
  # Adjust column widths
  width(j = "Gene Symbol", width = 1.2) %>%
  width(j = "Genomic Position", width = 1.5) %>%
  width(j = "Variant Type", width = 1.0) %>%
  width(j = "Allele Change", width = 1.0) %>%
  width(j = "Frequency (%)", width = 1.0) %>%
  width(j = "Mean VAF", width = 1.0) %>%
  width(j = "Population AF", width = 1.2) %>%
  width(j = "dbSNP ID", width = 1.5) %>%
  # Center align numeric columns
  align(j = c("Frequency (%)", "Mean VAF", "Population AF"), align = "center") %>%
  # Add table caption
  add_header_lines("Table 1. Recurrent Silent Mutations in HTLV-1 Cohort (>10% frequency)") %>%
  # Format font
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 11, part = "header") %>%
  # Add borders
  border_outer(border = fp_border(color = "black", width = 1)) %>%
  border_inner_h(border = fp_border(color = "gray", width = 0.5)) %>%
  border_inner_v(border = fp_border(color = "gray", width = 0.5))

# Save as Word document following project conventions
silent_genes_flextable %>%
  save_as_docx(path = "output/2508_silent_genes_publication_table.docx")

# Save as PNG for figures
silent_genes_flextable %>%
  save_as_image(path = "output/2508_silent_genes_publication_table.png", 
                width = 12, height = 8, res = 300)

# Create summary table for manuscript text
summary_stats_table <- silentgenes_csv %>%
  dplyr::summarise(
    total_variants = n(),
    genes_affected = n_distinct(Hugo_Symbol),
    mean_frequency = round(mean(Frequency), 1),
    median_frequency = round(median(Frequency), 1),
    frequency_range = paste0(round(min(Frequency), 1), "-", round(max(Frequency), 1)),
    mean_vaf = round(mean(mean_VAF), 3),
    variants_with_dbsnp = sum(!is.na(dbSNP_ID) & dbSNP_ID != ""),
    dbsnp_rate = round(variants_with_dbsnp / total_variants * 100, 1)
  )

# Export summary for manuscript
readr::write_csv(publication_table, "output/2508_silent_genes_publication_data.csv")
readr::write_csv(summary_stats_table, "output/2508_silent_genes_summary_stats.csv")

# Display the formatted table
silent_genes_flextable |>
  flextable::save_as_docx(path = "output/2508_silent_genes_publication_table.docx")


silentgenes_csv = silent_db2 |>
  filter(Hugo_Symbol!="Unknown") |>
  group_by(Hugo_Symbol, Chromosome, Start_Position, Variant_Type, POPAF, dbSNP_ID, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2) |>
  summarise(Frequency = n()/67*100, .groups = "drop",
            mean_VAF = mean(tumor_f, na.rm = TRUE)) |>
  filter(Frequency>10) |>
  arrange(desc(Frequency))
silentgenes_csv

# Create additional supplementary table with detailed annotations
supplementary_table <- silentgenes_csv %>%
  dplyr::left_join(
    silent_db2 %>% 
      dplyr::select(Hugo_Symbol, Start_Position, Variant_Classification, 
                   GO_Biological_Process, HGNC_Approved_name) %>%
      dplyr::distinct(),
    by = c("Hugo_Symbol", "Start_Position")
  ) %>%
  dplyr::select(
    Hugo_Symbol, HGNC_Approved_name, Start_Position, Variant_Classification,
    Frequency, mean_VAF, POPAF, dbSNP_ID, GO_Biological_Process
  ) %>%
  dplyr::rename(
    "Gene Symbol" = Hugo_Symbol,
    "Gene Name" = HGNC_Approved_name,
    "Genomic Position" = Start_Position,
    "Variant Classification" = Variant_Classification,
    "Frequency (%)" = Frequency,
    "Mean VAF" = mean_VAF,
    "Population AF" = POPAF,
    "dbSNP ID" = dbSNP_ID,
    "GO Biological Process" = GO_Biological_Process
  )

# Save supplementary table
supplementary_table %>%
  flextable() %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  save_as_docx(path = "output/2508_silent_genes_supplementary_table.docx")

cat("=== PUBLICATION TABLE SUMMARY ===\n")
print(summary_stats_table)

cat("\n=== TOP 10 RECURRENT SILENT MUTATIONS ===\n")
print(publication_table %>% head(15))


#---4- Correlation analysis betweenthe number of mutations and PVL
correlation_silentdb2 = silent_db2 |>
    mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
    mutate(ID=str_remove(ID, "^[0-9]+_")) |>
    group_by(ID) |>
    summarise(Frequency = n(), .groups = "drop") |>
    arrange(desc(Frequency)) |>
    right_join(metadata |> dplyr::select(PVLlog, Mean_HE, ID2, disease_2), by = c("ID"="ID2")) 
correlation_silentdb2
label = paste("Spearman rho =", 
                          round(cor(correlation_silentdb2$Frequency, correlation_silentdb2$PVLlog, 
                                    method="spearman", use="complete.obs"), 2),
                          "\np =", 
                          round(cor.test(correlation_silentdb2$Frequency, correlation_silentdb2$PVLlog, 
                                        method="spearman", exact = FALSE)$p.value, 3))


########### Mutation analysis
dataset2@data= dataset2@data |>
    dplyr::filter(!str_detect(Tumor_Sample_Barcode, "SS14|SS20|SS02|SS04|24_1339")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", "")) 

dataset2@variant.classification.summary = dataset2@variant.classification.summary |>
 dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

dataset2@variant.type.summary = dataset2@variant.type.summary |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

dataset2@clinical.data = dataset2@clinical.data |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_S\\d{1,3}_L\\d{3,4}__.*", "")) |>
    dplyr::mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "^\\d{2}_", ""))  

#### Anotate Clinical data

dataset2@clinical.data = dataset2@clinical.data|>
  left_join(metadata |> dplyr::select(ID2, PVLlog, Mean_HE, disease_2), 
            by = c("Tumor_Sample_Barcode"="ID2")) 

#### Test it
getSampleSummary(dataset2) |> View()
getClinicalData(dataset2)
library(grid)


("2508_maf_barplot.png", width = 2400, height = 1600, res = 300)
mafbarplot(dataset2, n=30, fontSize=0.7, legendfontSize = 0.7)
dev.off()

png("2508_maf_barplot2_summary.png", width = 2400, height = 1600, res = 300)
plotmafSummary(maf = dataset2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, textSize = c(1,1))
dev.off()

png("2508_sigpw.png", width = 2400, height = 1600, res = 300)
oncoplot(maf = dataset2, pathways = "sigpw", gene_mar = 8, fontSize = 0.7, topPathways = 3, collapsePathway = FALSE)
dev.off()

png("2508_smgpb.png", width = 2900, height = 1600, res = 300)
oncoplot(maf = dataset2, pathways = "smgbp", gene_mar = 8, fontSize = 0.7, topPathways = 3, collapsePathway = FALSE)
dev.off()
dataset2@clinical.data = dataset2@clinical.data |>
  dplyr::rename("log(PVL)" = PVLlog, "Mean HE" = Mean_HE, "Outcome disease" = disease_2)

dataset2@data |>
    dplyr::filter(!is.na(tumor_f)) |>
    dplyr::group_by(Hugo_Symbol) |>
    dplyr::summarise('mean VAF' = mean(tumor_f, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(desc('mean VAF')) -> genes_vaf
genes_vaf
### Include clinical data in oncoplot
tiff("2508_oncoplot_clinical2.tiff", width = 8, height = 5.2, units="in", res = 600)
oncoplot(maf = dataset2, 
         clinicalFeatures = c("Outcome disease", "log(PVL)", "Mean_HE"),
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         fontSize = 0.8,
         leftBarVlineCol = "blue",
         anno_height = 2,
         leftBarData = genes_vaf,
         leftBarLims = c(0, 0.1),
         legend_height=5)
dev.off()


#lollipop plot for KMT2DA, which is one of the most frequent mutated gene in Leukemia.
png("2506_KMT2DA.png", width = 2400, height = 1600, res = 200)
lollipopPlot(
  maf = dataset2,
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
  maf = dataset2,
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
  maf = dataset2,
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


### Plot Variant Allele Frequency (VAF) for the top 20 variants
png("2508_vaf.png", width = 2400, height = 1600, res = 300)
plotVaf(maf = dataset2, vafCol = "tumor_f", top = 20, )
dev.off()

### Plot Somatic Interactions
png("2508_Somatic.png", width = 2400, height = 1600, res = 300)
som = somaticInteractions(maf = dataset2, top = 25, pvalue = c(0.05, 0.01), fontSize=.6)
dev.off()


### Detecting cancer driver genes
mut = oncodrive(maf = dataset, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
mut
png("2506_Oncodrivers.png", width = 2400, height = 1600, res = 300)

drivers = mut |>
  select(Hugo_Symbol, pval, fdr, fract_muts_in_clusters) |>
  filter(fdr < 0.00017857, pval < 0.01, fract_muts_in_clusters >0.3) |>
  as.data.frame() #|>
  
drivers
  
  tidyplots::tidyplot(x=fract_muts_in_clusters, y=fdr) |>
  tidyplots::add_data_points_jitter()
plotOncodrive(res = mut, fdrCutOff = 0.05, useFraction = TRUE, labelSize = .4)
dev.off()
?plotOncodrive

#### Clinical enrichment analysis
fab.ce = clinicalEnrichment(maf = dataset, clinicalFeature = 'Disease')
fab.ce$groupwise_comparision[p_value < 0.05]


#### Heterogeneity
pacman::p_load(mclust)
het_h = inferHeterogeneity(maf = sub_maf, vafCol ="VAF")
het_h


#Requires BSgenome object
pak::pkg_install("BSgenome.Hsapiens.UCSC.hg38")

tnmatrix = trinucleotideMatrix(maf = dataset2, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
tnmatrix
plotApobecDiff(tnm = tnmatrix, maf = dataset2, pVal = 0.05)

pacman::p_load(NMF)
tnsignature = estimateSignatures(mat = tnmatrix, nTry = 6)
plotx = plotCophenetic(res = tnsignature)
plotx
maftools::plotSignatures(nmfRes = tnsignature, contributions = FALSE,)
tnsignature
extractSignatures(mat = tnsignature, n = 1)
