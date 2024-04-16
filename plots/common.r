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
eval_dir <- file.path(proj_dir, 'summaries/327/eval')
plots_dir <- file.path(proj_dir, 'plots')

cmrg_loci <- readLines(file.path(proj_dir, 'targets/327/CMRG_loci.txt'))

rgb <- function(r, g, b) {
    grDevices::rgb(r, g, b, maxColorValue = 255)
}

ascii_qvs <- function(x) {
     sub('–', '-', x) %>% sub('≥', '>=', .)
}

group_qvs <- function(summary, haps = T) {
    thresholds <- c(0, 17, 23, 33)
    labels <- c('<17', '17–23', '23–33', '≥33')
    get_qv_group <- function(qv) {
        n <- length(thresholds)
        findInterval(qv, thresholds)
    }
    
    filter(summary, !is.na(qv) & (query_type == 'hap') == haps) |>
    mutate(
        qv_cat = factor(labels[get_qv_group(qv)], levels = rev(labels)),
        qv_num = as.numeric(qv_cat),
    )
}

assign_filters <- function(summary) {
    cols <- list(warnings = '*', weight_dist = NA_real_,
        unexpl_reads = NA_real_, total_reads = 0.0)
    summary <- add_column(summary, !!!cols[setdiff(names(cols), colnames(summary))]) |>
    mutate(
        weight_dist = replace_na(weight_dist, 0),
        unexpl_reads = replace_na(unexpl_reads, 0),
        filt = warnings == '*'
            & weight_dist <= 30
            & (unexpl_reads <= 1000 | unexpl_reads <= total_reads * 0.2),
        filt = factor(ifelse(filt, 'Pass', 'Fail'), levels = c('Pass', 'Fail')),
    )
}

get_counts <- function(summary) {
    totals <- table(summary$locus)
    count(summary, locus, qv_cat, filt, name = 'count') |>
        complete(locus, qv_cat, filt, fill = list(count = 0)) |>
        mutate(
            locus = factor(locus),
            frac = as.numeric(count / totals[locus]),
            inter = interaction(qv_cat, filt, sep = '\n'),
            qv_num = as.numeric(qv_cat)
        )
}

bound_qv <- function(df, min_val = 40, edit0 = 0.5) {
    mutate(df, qv = ifelse(is.finite(qv), qv, pmax(min_val, -10 * log10(edit0 / size))))
}

.IMPORTED_COMMON <- T
