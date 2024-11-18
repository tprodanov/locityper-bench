library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringi)
library(crayon)
library(cowplot)
library(gtools)
library(Cairo)
cowplot::set_null_device('agg')

proj_dir <- '~/Data/proj/locityper'
eval_dir <- file.path(proj_dir, 'summaries/327/eval')
plots_dir <- file.path(proj_dir, 'plots')

cmrg_loci <- readLines(file.path(proj_dir, 'targets/327/CMRG_loci.txt'))

rgb <- function(r, g, b) {
    grDevices::rgb(r, g, b, maxColorValue = 255)
}

ascii_qvs <- function(x) {
    ifelse(x == 'No call', x,
        sub('–', '-', x) %>% sub('≥', '>=', .) %>% paste0('QV ', .)
    )
}

group_qvs <- function(summary, thresholds, haps = T) {
    stopifnot(thresholds[1] == 0)
    n <- length(thresholds)
    labels <- c(
        sprintf('<%d', thresholds[2]),
        sprintf('%d–%d', thresholds[2:(n-1)], thresholds[3:n]),
        sprintf('≥%d', thresholds[n]))
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
        qv_cat0 = qv_cat,
        qv_cat = ifelse(filt == 'Pass', as.character(qv_cat0), 'No call') |>
            factor(levels = c(levels(qv_cat0), 'No call')),
    )
}

get_counts <- function(summary) {
    totals <- table(summary$locus)
    count(summary, locus, qv_cat, name = 'count') |>
        complete(locus, qv_cat, fill = list(count = 0)) |>
        mutate(
            locus = factor(locus),
            frac = as.numeric(count / totals[locus]),
            qv_num = as.numeric(qv_cat)
        )
}

bound_qv <- function(df, min_val = 40, edit0 = 0.5, prefix = '') {
    #mutate(df, qv = ifelse(is.finite(qv), qv, pmax(min_val, -10 * log10(edit0 / size))))
    name_qv <- paste0(prefix, 'qv')
    name_size <- paste0(prefix, 'size')
    stopifnot(name_qv %in% colnames(df))
    stopifnot(name_size %in% colnames(df))
    df[[name_qv]] <- ifelse(is.finite(df[[name_qv]]),
        df[[name_qv]],
        pmax(min_val, -10 * log10(edit0 / df[[name_size]]))
        )
    df
}

draw_pdf <- function(filename, g, width = 15, height = 10, scale = 1, ...) {
    Cairo(type = 'pdf', file = filename,
        units = 'cm', width = width / scale, height = height / scale, ...)
    plot(g)
    while (!is.null(dev.list())) { dev.off() }
}

.IMPORTED_COMMON <- T
