library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringi)
library(crayon)
library(cowplot)
library(gtools)
cowplot::set_null_device('agg')

proj_dir <- '~/Data/proj/locityper'
plots_dir <- file.path(proj_dir, 'plots')

cmrg_loci <- readLines(file.path(proj_dir, 'targets/327/CMRG_loci.txt'))

rgb <- function(r, g, b) {
    grDevices::rgb(r, g, b, maxColorValue = 255)
}

preproc_summary <- function(summary, haps = T) {
    thresholds <- c(0, 17, 23, 33)
    labels <- c('<17', '17–23', '23–33', '≥33')
    get_qv_group <- function(qv) {
        n <- length(thresholds)
        findInterval(qv, thresholds)
    }

    summary <- filter(summary, (query_type == 'gt') != haps & !is.na(qv))
    cols <- list(warnings = '*', weight_dist = NA_real_,
        unexpl_reads = NA_real_, total_reads = 0.0)
    summary <- add_column(summary, !!!cols[setdiff(names(cols), colnames(summary))])

    mutate(summary,
        weight_dist = replace_na(weight_dist, 0),
        unexpl_reads = replace_na(unexpl_reads, 0),
        filt = warnings == '*'
            & weight_dist <= 30
            & (unexpl_reads < 1000 | unexpl_reads <= total_reads * 0.2),
        filt = factor(ifelse(filt, 'Pass', 'Fail'), levels = c('Pass', 'Fail')),
        qv_num = get_qv_group(qv),
        qv_cat = factor(labels[qv_num], levels = rev(labels))
    )
}

.IMPORTED_COMMON <- T
