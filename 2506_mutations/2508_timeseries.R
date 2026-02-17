#---1-Loading packages
pacman::p_load(maftools, tidyverse, data.table, tidyplots)
pak::pkg_install("corr")
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

dataset = merge_maf_objects('/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/Project_HTLV1 Peru_2025/2507_HTLVpaper/maf_htlv/series') 

totalmut = maftools::getSampleSummary(dataset) |>
    mutate(ID=str_remove(Tumor_Sample_Barcode, "_S[0-9]+.*$")) |>
    mutate(ID=str_remove(ID, "^[0-9]+_")) 

## Updating sheet4 with metadata and total tumor mutations

metadata = googlesheets4::read_sheet('https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=148794326#gid=148794326', sheet="Sheet4") 

metadata = metadata |>
    left_join(totalmut |> select(-Tumor_Sample_Barcode), by=c("ID1"="ID")) 

googlesheets4::sheet_write(metadata, ss='https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=148794326#gid=148794326', sheet="Sheet4")

### Re-reading the dataset
seriesdb = googlesheets4::read_sheet('https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=148794326#gid=148794326', sheet="Sheet4")

series1 = seriesdb |>
    select(ID2, PVL, time, Gini, total ) |>
    rename("PVL(%)"=PVL, "Oligoclonality Index (OCI)"=Gini, "Total Number of Mutations"=total) |>
    pivot_longer(-c(ID2, time), values_to = "score", names_to = "group") |>
    tidyplot(x = time, y = score, color=time) |> 
  add_line(group = ID2, color = "grey") |> 
  tidyplots::add_data_points(size=2) |>
  tidyplots::remove_y_axis_title() |>
  tidyplots::remove_x_axis_title() |>
  tidyplots::adjust_x_axis(labels=c("Baseline", "Follow-up")) |>
    tidyplots::remove_legend() |>
  add_test_pvalue() |>
    tidyplots::remove_caption() |>
    tidyplots::adjust_font(face="bold" ) |>
  tidyplots::split_plot(by=group) #|>
      tidyplots::save_plot("output/2508_timeseries.png", bg="transparent")

ggsave("2510_timeseries.tiff", series1, width = 7, height = 3, dpi = 600)


seriesdb |>
  select(time, total ) |>
  group_by(time) |>
  summarise(mean_total = mean(total, na.rm = TRUE))

matching_barcodes <- unique(dataset@data$Tumor_Sample_Barcode[str_detect(dataset@data$Tumor_Sample_Barcode, "IRID015|IRID039|IRID062|IRID089")])
matching_barcodes2 = unique(dataset@data$Tumor_Sample_Barcode[str_detect(dataset@data$Tumor_Sample_Barcode, "SS04|SS14|SS20|SS02")])

sub1 = maftools::subsetMaf(dataset, tsb = matching_barcodes)
sub2 = maftools::subsetMaf(dataset, tsb = matching_barcodes2)
tiff("2510_timeseries_maf_barplot.tiff", width = 6, height = 4, res = 600, unit="in")
coBarplot(m1 = sub2, m2 = sub1, m1Name = "Baseline", m2Name = "Follow-up", geneSize=.6, legendTxtSize = .7, showPct = FALSE, titleSize = 1)
dev.off()
