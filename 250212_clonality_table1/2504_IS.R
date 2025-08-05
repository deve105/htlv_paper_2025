## Packages
pacman::p_load(
    tidyverse, tidyplots
)

## Loading file
overall_coverage = read_tsv("./data/2504_mosdepth_global.tsv", col_names=FALSE)

overall_cov = overall_coverage |>
    dplyr::rename(chrom=1, length=2, bases=3, median=4, ID=7) |>
    dplyr::mutate(median=as.numeric(median), 
        less=ifelse(as.numeric(median)<5, "LC", "OK"),
        peru=ifelse(str_detect(ID,"IRID|SS"), "Peru", "Japan"),
        ID=str_replace_all(ID, "2._", "")) |>
        dplyr::filter(chrom=="total", peru=="Peru") 

class(overall_cov$median)
overall_cov |>
    dplyr::filter(median<5)


overall_cov |>
    summarise(
    n = n(),
    mean = mean(median, na.rm = TRUE),
    sd = sd(median, na.rm = TRUE),
    se = sd / sqrt(n),
    lower_ci = mean - qt(0.975, df = n - 1) * se,  # 95% CI lower bound
    upper_ci = mean + qt(0.975, df = n - 1) * se   # 95% CI upper bound
  )


overall_cov |> 
    summarise(
        n = n(),
        medianx = median(median, na.rm = TRUE),
        Q1 = quantile(median, 0.1, na.rm = TRUE),
        Q3 = quantile(median, 0.7, na.rm = TRUE),
        IQR = IQR(median, na.rm = TRUE)
    )

overall_cov |>
        dplyr::mutate(median=median+1 ) |>
        tidyplots::tidyplot(y=ID, x=median) |>
        tidyplots::add_median_bar(width=.5) |>
        tidyplots::sort_y_axis_labels(.reverse=TRUE) |>
        tidyplots::adjust_y_axis_title("Sample ID", face="bold") |>
        tidyplots::adjust_x_axis_title("Whole HTLV-1 Coverage per Sample\n log10(Coverage + 1)", face="bold") |>
        tidyplots::adjust_x_axis(transform = "log10", labels = scales::trans_format("log10", scales::math_format(10^.x))) |>
        tidyplots::adjust_font(fontsize = 6) |>
        tidyplots::add_reference_lines(x = 159, linetype = "dashed", linewidth = .5, color="red") |>
        tidyplots::adjust_size(width = 5, height = 13, unit = "cm") |>
        tidyplots::save_plot("2504_HTLV1coverage.png", bg="transparent")
