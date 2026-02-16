
pacman::p_load(tidyverse)

# Load the VCF

htlvis = read.delim("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/2507_concatenatedtab_vcf.tsv", header = TRUE, sep = "\t") |>
    dplyr::filter(str_detect(X.CHROM, "adjusted")) |>
    mutate(VAF = str_extract(htlvis$"X22_SS01_S238_L008__sorted.bam.viral.bam", ":[^:]*:") |>
    str_remove_all(":") |>
    as.numeric()) |>
    rename(sample=X22_SS01)

colnames(htlvis)

htlvis$"X22_SS01_S238_L008__sorted.bam.viral.bam"
str_extract(htlvis$"X22_SS01_S238_L008__sorted.bam.viral.bam", ":[^:]*:") |>
    str_remove_all(":") |>
    as.numeric()

htlvis$POS
plot(htlvis$POS)

pacman::p_load(DescTools)

# For every sample, calculate Gini index from VAF
gini_by_sample <- htlvis %>%
  group_by(sample) %>%
  summarise(gini_index = Gini(VAF, na.rm = TRUE))
DT::datatable(gini_by_sample)
