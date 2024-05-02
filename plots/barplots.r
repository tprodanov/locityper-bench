if (!exists('.IMPORTED_COMMON')) { source('common.r') }

draw_barplot <- function(counts, panel_width = 71, colors = NULL,
    fill_label = 'Accuracy (QV)')
{
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
        scale_y_continuous('Fraction of haplotypes   ', expand = c(0, 0)) +
        scale_fill_manual(paste0(fill_label, ' '),
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

aggr_summaries <- function(summaries, key) {
    total <- nrow(summaries)
    aggrs <- mutate(summaries,
        qv_cat = factor(ifelse(filt == 'Fail', 'N/A', as.character(qv_cat)),
            levels = c(levels(summaries$qv_cat), 'N/A'))) |>
        count(qv_cat, name = 'count') |>
        mutate(key = key, frac = count / total)
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
            ascii_qvs(levels(counts$qv_cat)[i]),
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
    if (total_fail != 0 & total_pass != 0) {
        subdescribe(counts_pass, F)
    }
    cat(sprintf('Failed %s\n', percentage(total_fail, total)))
    if (total_fail != 0 & total_pass != 0) {
        subdescribe(counts_fail, F)
    }
}
describe_several <- function(counts, pattern = '') {
    for (key in names(counts)) {
        if (grepl(pattern, key)) {
            cat(sprintf('\n=== %s ===\n', bold(key)))
            describe(counts[[key]])
        }
    }
}

#####################################

ver <- 'v0.14.2-4'
summaries <- list()
counts <- list()
for (tech in c('illumina', 'hifi', 'ont', 'sim')) {
    for (loo in c(F, T)) {
        key <- sprintf('%s%s', tech, ifelse(loo, '_loo', ''))
        summaries[[key]] <- read.csv(
            sprintf('%s/%s_%s_%s.csv.gz',
                eval_dir, ifelse(loo, 'loo', 'full'), tech, ver),
            sep = '\t', comment = '#') |>
            filter(locus %in% cmrg_loci) |>
            group_qvs() |>
            bound_qv() |>
            assign_filters()
        counts[[key]] <- get_counts(summaries[[key]])
    }
}

summaries$nygc <- read.csv(
    sprintf('%s/comparison/NYGC/summaries/nygc.eval.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs() |> assign_filters()
counts$nygc <- get_counts(summaries$nygc)
describe(counts$nygc)

summaries$avail <- mutate(summaries$illumina_loo,
    edit = avail_edit,
    size = avail_size,
    div = edit / size,
    qv = avail_qv, filt = 'Pass') |>
    group_qvs() |>
    bound_qv()
counts$avail <- get_counts(summaries$avail)
describe(counts$avail)

local({
    summary_1kgp <- read.csv(file.path(eval_dir, '../full_1kgp_v0.14.4-0.csv.gz'),
        sep = '\t', comment = '#') |>
        filter(locus %in% cmrg_loci) |>
        assign_filters() |>
        select(sample, locus, filt)
    summaries$trios <<- read.csv(file.path(eval_dir, 'trio_conc.csv.gz'),
            sep = '\t', comment = '#') |>
        filter(locus %in% cmrg_loci) |>
        group_qvs(haps = T) |>
        bound_qv() |>
        left_join(summary_1kgp, by = join_by(indiv == sample, locus)) |>
        left_join(summary_1kgp, by = join_by(mother == sample, locus),
            suffix = c('', '.mother')) |>
        left_join(summary_1kgp, by = join_by(father == sample, locus),
            suffix = c('', '.father')) |>
        mutate(
            filt = ifelse(filt == 'Pass' & filt.mother == 'Pass' & filt.father == 'Pass',
                'Pass', 'Fail') |> factor(levels = c('Pass', 'Fail'))
        )
    summaries$trios2 <<- mutate(summaries$trios, filt = 'Pass')
    counts$trios <<- get_counts(summaries$trios)
    counts$trios2 <<- get_counts(summaries$trios2)
})

barplots <- list()
for (key in names(counts)) {
    barplots[[key]] <- draw_barplot(counts[[key]])
}

##########################

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
(lost_qv_panel <- ggplot(lost_qvs) +
    geom_blank() +
    geom_bar(aes(lost_qv + 0.5, cum_frac, fill = 'Cumulative'),
        alpha = 0.4,
        stat = 'identity', width = 1) +
    geom_bar(aes(lost_qv + 0.5, frac, fill = 'Histogram'),
        stat = 'identity', width = 1) +
    facet_wrap(~ tech, ncol = 2) +
    coord_cartesian(xlim = c(0, 10)) +
    scale_x_continuous('Lost accuracy (QV)', expand = expansion(mult = 0.00),
        breaks = c(0, 5, 10)) +
    scale_y_continuous('Fraction of haplotypes', expand = expansion(mult = 0.00),
        breaks = seq(0, 1, 0.2), minor_breaks = c(0.9)) +
    scale_fill_manual(NULL,
        breaks = c('Histogram', 'Cumulative'),
        values = c(rgb(59, 61, 126), '#1961AC00')) +
    guides(fill = guide_legend(override.aes = list(color = '#000000C0', linewidth = 0.2))) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linewidth = 0.2, color = 'gray50', linetype = '22'),
        panel.grid.major = element_line(linewidth = 0.2, color = 'gray50', linetype = '22'),
        strip.background = element_rect(color = NA, fill = 'gray100'),
        strip.text = element_text(margin = margin(l = 4, t = 1, b = 1), hjust = 0),
        panel.spacing.x = unit(0.7, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -5, r = 6),
        legend.key.height = unit(0.9, 'lines'),
        legend.key.width = unit(0.4, 'lines'),
        legend.text = element_text(margin = margin(l = -3)),
        plot.margin = margin(l = 14, r = 6, t = 5, b = 2),
    ))

aggrs <- lapply(names(summaries),
    function(name) aggr_summaries(summaries[[name]], name)) %>%
    do.call(rbind, .) |>
    filter(key != 'trios') |>
    mutate(
        panel = case_when(
            grepl('^(illumina|hifi|ont|sim|nygc)$', key) ~ 'a',
            grepl('^(avail|.*_loo)$', key) ~ 'b',
            key == 'trios2' ~ 'c'),
        facet = case_when(
            grepl('^(illumina|hifi|ont|sim)$', key) ~ 'Locityper ●',
            grepl('^(illumina|hifi|ont|sim)_loo$', key) ~ 'Locityper ○',
            key == 'avail' ~ ' ',
            key == 'nygc' ~ ' ',
            key == 'trios2' ~ ' ',
            T ~ '── Locityper ──') |>
            factor(levels = c('Locityper ●', 'Locityper ○', ' ')),
        label = case_when(
            grepl('illumina', key) ~ 'Illumina',
            grepl('sim', key) ~ 'Simulated',
            grepl('hifi', key) ~ 'HiFi',
            grepl('ont', key) ~ 'ONT',
            key == 'avail' ~ 'Avail.',
            key == 'nygc' ~ '1KGP',
            key == 'trios2' ~ ' Trio conc.') |>
            factor(levels = c('Illumina', 'Simulated', 'HiFi', 'ONT',
                'Avail.', '1KGP', ' Trio conc.'))
    )

aggr_plots <- list()
for (letter in c('a', 'b', 'c')) {
    aggr_plots[[letter]] <- filter(aggrs, panel == letter) |>
        ggplot() +
        geom_bar(aes(label, frac, fill = qv_cat),
            position = position_stack(reverse = T),
            stat = 'identity', width = 1, color = '#000000C0', linewidth = 0.1) +
        facet_grid(. ~ facet, space = 'free', scales = 'free') +
        scale_fill_manual('     QV  ',
            breaks = c(rev(levels(summaries$illumina$qv_cat)), 'N/A'),
            values = c('#ec4900', '#e1a900', '#9ebe91', '#3b9ab2', 'gray80'),
            drop = F) +
        geom_hline(yintercept = seq(0.2, 0.8, 0.2), linetype = '22',
            color = 'black', alpha = 0.7, linewidth = 0.2) +
        geom_hline(yintercept = c(0, 1),
            color = 'black', alpha = 0.7, linewidth = 0.2) +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous('Fraction of haplotypes',
            breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            panel.background = element_rect(fill = NA),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'bottom',
            legend.margin = margin(t = 1, b = 5),
            legend.title = element_text(size = 11, family = 'Overpass SemiBold'),
            legend.text = element_text(family = 'Overpass',
                margin = margin(l = -2, r = 3), hjust = 0.5),
            legend.key.size = unit(0.9, 'lines'),
            strip.background = element_rect(color = NA, fill = 'white'),
            strip.text.x = element_text(margin = margin(t = 2, b = 2)),
            panel.spacing.x = unit(0.5, 'lines'),
            axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
            plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(l = 5, r = 5, t = 8, b = 7)
        )
}

left_side <- plot_grid(
    plot_grid(
        aggr_plots$a + theme(legend.position = 'none'),
        aggr_plots$b + theme(legend.position = 'none', axis.title.y = element_blank()),
        aggr_plots$c + theme(legend.position = 'none', axis.title.y = element_blank()),
        nrow = 1,
        rel_widths = c(0.45, 0.40, 0.15)
    ),
    get_legend(aggr_plots$a),
    nrow = 2,
    rel_heights = c(0.95, 0.05)
    )
plot_grid(
    left_side,
    lost_qv_panel,
    rel_widths = c(0.64, 0.36),
    nrow = 1)
ggsave(file.path(plots_dir, 'fig2.png'), bg = 'white',
    width = 8.1, height = 5, dpi = 600, scale = 0.7)

ggsave(file.path(plots_dir, 'fig2_left.png'), left_side, bg = 'white',
    width = 4.5, height = 6, dpi = 600, scale = 0.8)

#####################

ggsave(sprintf('%s/barplots/illumina.png', plots_dir),
    barplots$illumina,
    width = 10, height = 7, dpi = 600, scale = 0.75)
plot_grid(
    barplots$sim + labs(subtitle = 'Locityper, Simulated data'),
    barplots$nygc + labs(subtitle = '1KGP phased call set'),
    nrow = 2,
    labels = c('a', 'b'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )
ggsave(sprintf('%s/barplots/supp1ab.png', plots_dir),
        width = 10, height = 13, dpi = 600, scale = 0.75)

plot_grid(
    barplots$hifi + labs(subtitle = 'Locityper, PacBio HiFi'),
    barplots$ont + labs(subtitle = 'Locityper, Oxford Nanopore'),
    nrow = 2,
    labels = c('c', 'd'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )
ggsave(sprintf('%s/barplots/supp1cd.png', plots_dir),
    width = 10, height = 13, dpi = 600, scale = 0.75)

plot_grid(
    barplots$illumina_loo + labs(subtitle = 'Locityper, Illumina, leave-one-out'),
    barplots$sim_loo + labs(subtitle = 'Locityper, Simulated data, leave-one-out'),
    nrow = 2,
    labels = c('a', 'b'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )
ggsave(sprintf('%s/barplots/supp2ab.png', plots_dir),
    width = 10, height = 13, dpi = 600, scale = 0.75)
plot_grid(
    barplots$hifi_loo + labs(subtitle = 'Locityper, PacBio HiFi, leave-one-out'),
    barplots$ont_loo + labs(subtitle = 'Locityper, Oxford Nanopore, leave-one-out'),
    nrow = 2,
    labels = c('c', 'd'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )
ggsave(sprintf('%s/barplots/supp2cd.png', plots_dir),
    width = 10, height = 13, dpi = 600, scale = 0.75)

ggsave(sprintf('%s/barplots/avail.png', plots_dir),
    draw_barplot(counts$avail, fill_label = 'Availability (QV)'),
    width = 10, height = 6.5, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/trios.png', plots_dir),
    draw_barplot(counts$trios, fill_label = 'Concordance (QV)'),
    width = 10, height = 6.5, dpi = 600, scale = 0.75)

##################

filter(summaries$illumina) |>
    group_by(locus) |>
    filter(mean(qv_num <= 2) >= 0.95) |>
    with(length(unique(locus)))

describe_several(counts, 'hifi_loo|avail')
describe_several(counts, 'ont')

filter(lost_qvs, lost_qv == 4 | lost_qv == 9) |> select(!c(n, frac))
filter(summaries$illumina_loo, filt == 'Pass') |>
    mutate(lost_qv = pmin(avail_qv, 33) - pmin(qv, 33)) |>
    with(mean(lost_qv < 5))

with(summaries$trios, mean(qv))
aggregate(qv ~ locus, summaries$trios, mean) |> with(min(qv))

trios_union <- union(summaries$trios$indiv, summaries$trios$mother) |>
    union(summaries$trios$father)

hprc_samples <- readLines('~/Data/data/samples/HPRC.txt')
setdiff(hprc_samples, trios_union) |> paste0(collapse = ', ')
describe(counts$trios)
