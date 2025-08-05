####-----------------------------------------------
pacman::p_load(maftools, tidyverse, data.table)

####-----------------------------------------------
#### Merging MAF files from Mutect2 and Funcotator
#---1-Loading Packages
pacman::p_load(maftools, tidyverse, data.table)

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
###

#---metadata oficial
metadata = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=1412798579#gid=1412798579", sheet="Sheet3") |>
    mutate(status1 = ifelse(Status == "vivo", "alive", "dead/lost" )) |>
    mutate(PVLlog = as.numeric(log(PVL)))
#---MAF official
# This code merges MAF files from a specified directory, ensuring that each file's sample name
dataset = merge_maf_objects("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/2507_HTLVpaper/maf_htlv") 

## Correlation analysis betweenthe number of mutations and PVL
# This code calculates the correlation between the total number of mutations and the PVL log value,
# and visualizes it using a scatter plot with a linear regression line.
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

library(tidyplots)
correlation |>
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
    tidyplots::adjust_font(face="bold") |>
    tidyplots::save_plot("output/2508_PVL_vs_total_mutations.png", bg="transparent")

## MAF summary
png("2508_summary_maf_barplot.png", width = 1800, height = 2400, res = 300)
plotmafSummary(maf = dataset, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, top = 20, titvRaw = FALSE)
dev.off()

getSampleSummary(dataset) |>
  dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) |>
  t() |>



col_sums
    

colnames(dataset@data)[which(str_detect(colnames(dataset@data), "ClinVar"))]

dataset@data$ClinVar_VCF_ID

  # Combine all data tables
  merged_data <- rbindlist(all_maf_data, use.names = TRUE, fill = TRUE)
  # Create final merged MAF object
  # (Message suppressed)
  merged_maf <- read.maf(maf = merged_data, verbose = FALSE)
  # (Summary statistics and print messages suppressed)
  # Save merged MAF object if output file specified
  if (!is.null(output_file)) {
    write.mafSummary(maf = merged_maf, basename = tools::file_path_sans_ext(output_file))
    fwrite(merged_maf@data, output_file, sep = "\t", quote = FALSE, na = "")
    # cat("\nMerged MAF file written to:", output_file, "\n")
  }
  return(merged_maf)
}

## Silent mutations calculation
silent_maf_objects <- function(maf_directory = ".", output_file = NULL) {
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
      maf_data <- maf_obj@maf.silent
      maf_data$Tumor_Sample_Barcode <- sample_name
      maf_objects[[sample_name]] <- maf_data
    }, error = function(e) {
      warning(paste("Error reading", sample_name, ":", e$message))
    })
  }
  if (length(maf_objects) == 0) {
    stop("No valid MAF files could be processed")
  }
  merged_data <- rbindlist(maf_objects, use.names = TRUE, fill = TRUE)
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
