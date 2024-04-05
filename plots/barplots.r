if (!exists('.IMPORTED_COMMON')) { source('common.r') }

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

draw_barplot <- function(counts, panel_width = 71, colors = NULL) {
    n_loci <- length(levels(counts$locus))
    panels <- ceiling(n_loci / panel_width)
    panel_size <- ceiling(n_loci / panels)
    counts$panel <- (as.numeric(counts$locus) - 1) %/% panel_size

    any_lowqual <- any(counts$filt == 'Fail' & counts$count > 0)
    if (is.null(colors)) {
        # From wesanderson::wes_palette("Zissou1", N, type = "continuous") with N = 10 and 51.
        colors <- c('#3b9ab2', '#9ebe91', '#e1a900', '#ec4900')
    }

    ggplot(counts) +
        geom_bar(aes(locus, frac, group = inter, fill = qv_cat, alpha = filt),
            stat = 'identity', position = position_stack(reverse = T),
            width = 1, color = '#000000C0', linewidth = 0.1) +
        facet_wrap(~ panel, ncol = 1, scales = 'free_x') +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous('Fraction of samples   ', expand = c(0, 0)) +
        scale_fill_manual('Accuracy (QV) ',
            breaks = rev(levels(counts$qv_cat)),
            values = rev(colors)
            ) +
        scale_alpha_manual(' Filter ',
            breaks = c('Fail', 'Pass'),
            values = c(0.3, 1)) +
        guides(
            fill = guide_legend(order = 1),
            alpha = if (any_lowqual) { guide_legend(order = 2) } else { 'none' }) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(linetype = '22',
                color = '#000000C0', linewidth = 0.2),
            panel.background = element_rect(fill = NA),
            panel.ontop = TRUE,
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            legend.position = 'bottom',
            legend.margin = margin(t = -10),
            legend.title = element_text(size = 11, family = 'Overpass SemiBold'),
            legend.text = element_text(family = 'Overpass', margin = margin(l = -2, r = 3), hjust = 0.5),
            legend.key.size = unit(0.9, 'lines'),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing.y = unit(0.1, 'lines'),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(
                size = 6.2, face = 'italic', angle = 45, vjust = 1, hjust = 1)
        )
}

percentage <- function(a, b) {
    sprintf('%5d/%5d (%5.1f%%)', a, b, 100 * a / b)
}
subdescribe <- function(counts, header) {
    n <- sum(counts$count)
    if (header) {
        cat(sprintf('%-12s  %-20s  %-20s  %-20s\n',
            '', 'this accuracy', 'this or more acc.', 'this or less acc.'))
    }
    for (i in 1:length(levels(counts$qv_cat))) {
        cat(sprintf('    QV %-5s  %20s  %20s  %20s\n',
            stri_enc_toascii(levels(counts$qv_cat)[i]),
            percentage(with(counts, sum((qv_num == i) * count)), n),
            percentage(with(counts, sum((qv_num <= i) * count)), n),
            percentage(with(counts, sum((qv_num >= i) * count)), n)
            ))
    }
}
describe <- function(counts) {
    subdescribe(counts, T)
    counts_pass <- filter(counts, filt == 'Pass')
    counts_fail <- filter(counts, filt == 'Fail')
    total <- sum(counts$count)
    total_pass <- sum(counts_pass$count)
    total_fail <- sum(counts_fail$count)

    cat(sprintf('Passed %s\n', percentage(total_pass, total)))
    if (total_fail * total_pass != 0) {
        subdescribe(counts_pass, F)
    }
    cat(sprintf('Failed %s\n', percentage(total_fail, total)))
    if (total_fail * total_pass != 0) {
        subdescribe(counts_fail, F)
    }
}

#####################################

eval_dir <- file.path(proj_dir, 'summaries/eval')
ver <- 'v0.14.2-4'

summaries <- list()
counts <- list()
for (tech in c('illumina', 'hifi', 'ont', 'sim')) {
    for (loo in c(F, T)) {
        key <- sprintf('%s%s', tech, ifelse(loo, '_loo', ''))
        summaries[[key]] <- read.csv(
            sprintf('%s/327%s_%s_%s.csv.gz',
                eval_dir, ifelse(loo, 'LOO', ''), tech, ver),
            sep = '\t', comment = '#') |>
            filter(locus %in% cmrg_loci) |>
            preproc_summary()
        counts[[key]] <- get_counts(summaries[[key]])
    }
}

for (key in names(counts)) {
    cat(sprintf('\n=== %s ===\n', bold(key)))
    describe(counts[[key]])
}
barplots <- list()
for (key in names(counts)) {
    barplots[[key]] <- draw_barplot(counts[[key]])
}

summaries$nygc <- read.csv(
    sprintf('%s/comparison/NYGC/summaries/nygc.eval.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    preproc_summary()
counts$nygc <- get_counts(summaries$nygc)
barplots$nygc <- draw_barplot(counts$nygc)
describe(counts$nygc)

##########################

summaries$avail <- mutate(summaries$hifi_loo, qv = avail_qv, filt = 'Pass')
counts$avail <- get_counts(summaries$avail)
{
    aggrs <- data.frame()
    for (key in names(summaries)) {
        total <- nrow(summaries[[key]])
        if (is.null(summaries[[key]]$avail_qv)) {
            summaries[[key]]$avail_qv <- Inf
        }
        aggrs <- mutate(summaries[[key]],
                qv_cat = case_when(
                    filt == 'Fail' ~ 'N/A',
                    T ~ as.character(qv_cat)
                )
            ) |>
            count(qv_cat, name = 'count') |>
            mutate(key = key, frac = count / total) %>%
            rbind(aggrs, .)
    }
    aggrs <- mutate(aggrs,
        qv_cat = factor(qv_cat,
            levels = c(levels(counts$illumina$qv_cat), 'N/A')),
        subfigure = ifelse(endsWith(key, '_loo') | key == 'avail', 'b', 'a'),
        panel = factor(case_when(
            key == 'nygc' ~ ' ',
            key == 'avail' ~ '   ',
            # endsWith(key, '_loo') ~ 'Locityper, LOO',
            T ~ '── Locityper ──'), levels = c('   ', '── Locityper ──', ' ')),
        tech = factor(case_when(
            startsWith(key, 'illumina') ~ 'Illumina',
            startsWith(key, 'sim') ~ 'Simulated',
            startsWith(key, 'hifi') ~ 'HiFi',
            startsWith(key, 'ont') ~ 'ONT',
            key == 'avail' ~ 'Available',
            key == 'nygc' ~ '1KGP'),
            levels = c('Illumina', 'Simulated', 'HiFi', 'ONT', 'Available', '1KGP'))
        )
}
aggr_plots <- list()
for (label in c('a', 'b')) {
    aggr_plots[[label]] <- filter(aggrs, subfigure == label) |>
        ggplot() +
        geom_bar(aes(tech, frac, fill = qv_cat),
            position = position_stack(reverse = T),
            stat = 'identity', width = 1, color = '#000000C0', linewidth = 0.1) +
        facet_grid(. ~ panel, space = 'free', scales = 'free') +
        scale_fill_manual('QV',
            breaks = c(rev(levels(counts$illumina$qv_cat)), 'N/A'),
            values = c('#ec4900', '#e1a900', '#9ebe91', '#3b9ab2', 'gray80'),
            drop = F) +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous('Fraction of haplotypes', expand = c(0, 0)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(linetype = '22',
                color = '#000000C0', linewidth = 0.2),
            panel.background = element_rect(fill = NA),
            panel.ontop = TRUE,
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            legend.position = 'right',
            legend.margin = margin(b = -7),
            legend.title = element_text(size = 11, family = 'Overpass SemiBold'),
            legend.text = element_text(family = 'Overpass',
                margin = margin(l = -2, r = 3), hjust = 0.5),
            legend.key.size = unit(0.9, 'lines'),
            strip.background = element_rect(color = NA, fill = 'white'),
            strip.text.x = element_text(margin = margin(t = 2, b = 2)),
            panel.spacing.x = unit(0.5, 'lines'),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
        )
}

{
    lost_qvs <- data.frame()
    MAX_QV <- 33
    for (tech in c('illumina', 'hifi', 'ont', 'sim')) {
        lost_qvs <- filter(summaries[[paste0(tech, '_loo')]], filt == 'Pass') |>
            mutate(lost_qv = floor(pmin(avail_qv, MAX_QV) - pmin(qv, MAX_QV))) |>
            count(lost_qv) |>
            complete(lost_qv = 0:MAX_QV, fill = list(n = 0)) |>
            mutate(frac = n / sum(n), cum_frac = cumsum(frac), tech = tech) %>%
            rbind(lost_qvs, .)
    }
    lost_qvs <- mutate(lost_qvs,
        tech = recode_factor(tech,
            'illumina' = 'Illumina', 'sim' = 'Simulated',
            'hifi' = 'HiFi', 'ont' = 'ONT'))
}
lost_qv_panel <- ggplot(lost_qvs) +
    geom_blank() +
    geom_bar(aes(lost_qv + 0.5, cum_frac),
        alpha = 0.5, fill = rgb(25, 97, 172),
        stat = 'identity', width = 1) +
    geom_bar(aes(lost_qv + 0.5, frac),
        fill = rgb(59, 61, 126),
        stat = 'identity', width = 1) +
    geom_text(aes(x = Inf, y = Inf, label = 'Cumulative \nfraction '),
        family = 'Overpass', size = 2.5, angle = 0, hjust = 1, vjust = 1.5,
        data = data.frame(tech = factor('Illumina', levels = levels(lost_qvs$tech)))) +
    facet_wrap(~ tech, ncol = 2) +
    coord_cartesian(xlim = c(0, 10)) +
    scale_x_continuous('Lost accuracy (QV)', expand = expansion(mult = 0.00),
        minor_breaks = NULL, breaks = c(0, 5, 10)) +
    scale_y_continuous('Fraction of haplotypes', expand = expansion(mult = 0.00),
        breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color = NA, fill = 'gray100'),
        strip.text = element_text(margin = margin(l = 5, t = 1, b = 1), hjust = 0),
        panel.spacing.x = unit(0.7, 'lines'),
    )

plot_grid(
    aggr_plots$a + theme(legend.position = 'none'),
    aggr_plots$b,
    lost_qv_panel,
    rel_widths = c(0.25, 0.35, 0.4),
    nrow = 1,
    labels = letters,
    label_fontfamily = 'Overpass',
    label_fontface = 'bold',
    label_y = 1.01
    )
ggsave('~/Downloads/1.png',
    width = 10, height = 5, dpi = 600, scale = 0.85)
