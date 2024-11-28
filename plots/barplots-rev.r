if (!exists('.IMPORTED_COMMON')) { source('common.r') }

zissou <- wesanderson::wes_palette("Zissou1", 20, type = "continuous")
zissou100 <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
scales::show_col(zissou100, F)
THRESHOLDS <- c(0, 17, 23, 33, 43)
COLORS <- zissou100[c(1, 11, 35, 65, 92)]
# scales::show_col(COLORS, F)
stopifnot(length(THRESHOLDS) == length(COLORS))

draw_barplot <- function(summary, panel_width = 71, fill_label = 'Accuracy (QV)')
{
    counts <- get_counts(summary)
    n_loci <- length(levels(counts$locus))
    panels <- ceiling(n_loci / panel_width)
    panel_size <- ceiling(n_loci / panels)
    counts$panel <- (as.numeric(counts$locus) - 1) %/% panel_size

    ggplot(counts) +
        geom_bar(aes(locus, frac, fill = qv_cat),
            stat = 'identity', position = position_stack(reverse = F),
            width = 1, color = '#000000C0', linewidth = 0.1) +
        facet_wrap(~ panel, ncol = 1, scales = 'free_x') +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous('Fraction of haplotypes   ', expand = c(0, 0)) +
        scale_fill_manual(paste0(fill_label, ' '),
            breaks = levels(counts$qv_cat),
            values = c(COLORS, 'gray80')) +
        guides(fill = guide_legend(nrow = 1)) +
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
            legend.text = element_text(family = 'Overpass',
                margin = margin(l = 2, r = 3)),
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
    count(summaries, qv_cat, name = 'count') |>
        mutate(key = key, frac = count / total)
}

percentage <- function(a, b) {
    sprintf('%5d/%5d (%5.1f%%)', a, b, 100 * a / b)
}
.subdescribe <- function(summary) {
    cat(sprintf('    Mean QV: %.1f,  Median: %.1f\n', mean(summary$qv), median(summary$qv)))
    cats <- levels(summary$qv_cat)
    n <- nrow(summary)
    x <- summary$qv_num
    for (i in 1:length(cats)) {
        if (cats[i] != 'No call') {
            cat(sprintf('    %-8s  %20s  %20s  %20s\n',
                ascii_qvs(cats[i]),
                percentage(sum(x == i), n),
                percentage(sum(x <= i), n),
                percentage(sum(x >= i), n)
                ))
        }
    }
}
describe <- function(summary, short = T) {
    if (!any(summary$qv_cat == 'No call')) {
        .subdescribe(summary)
    } else {
        summary1 <- filter(summary, qv_cat != 'No call')
        cat(sprintf('Passed filtering: %s\n', percentage(nrow(summary1), nrow(summary))))
        .subdescribe(summary1)
        if (!short) {
            summary2 <- filter(summary, qv_cat == 'No call') |>
                mutate(qv_cat = qv_cat0, qv_num = qv_num0)
            summary3 <- mutate(summary, qv_cat = qv_cat0, qv_num = qv_num0)
            cat(sprintf('No call: %s\n', percentage(nrow(summary2), nrow(summary))))
            .subdescribe(summary2)
            cat('All:\n')
            .subdescribe(summary3)
        }
    }
}
describe_several <- function(summaries, pattern = '') {
    for (key in names(summaries)) {
        if (grepl(pattern, key)) {
            cat(sprintf('\n=== %s ===\n', bold(key)))
            describe(summaries[[key]])
        }
    }
}

#####################################

s <- read.csv(sprintf('%s/loo_illumina_v0.17.0.csv.gz', eval_dir),
    sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs(THRESHOLDS) |>
    bound_qv()
assign_filters(s) |> describe(F)
describe(s)

#################

vers <- c('v0.17.0') #, 'v0.14.2-4')
summaries <- list()
for (tech in c('illumina', 'hifi', 'ont', 'sim')) {
    for (loo in c(F, T)) {
        for (ver in vers) {
            filename <- sprintf('%s/%s_%s_%s.csv.gz',
                eval_dir, ifelse(loo, 'loo', 'full'), tech, ver)
            if (!file.exists(filename)) {
                next
            }
            key <- sprintf('%s%s', tech, ifelse(loo, '_loo', ''))
            summaries[[key]] <- read.csv(
                sprintf('%s/%s_%s_%s.csv.gz',
                    eval_dir, ifelse(loo, 'loo', 'full'), tech, ver),
                sep = '\t', comment = '#') |>
                filter(locus %in% cmrg_loci) |>
                group_qvs(THRESHOLDS) |>
                bound_qv() |>
                assign_filters()
            break
        }
    }
}
summaries$nygc <- read.csv(
    sprintf('%s/comparison/NYGC/summaries/nygc.eval.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs(THRESHOLDS) |> bound_qv() |>
    mutate(filt = 'Pass') # |> assign_filters()
describe(summaries$nygc)

summaries$sniffles <- read.csv(
    sprintf('%s/hifi-phasing/summaries/sniffles.eval.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs(THRESHOLDS) |> bound_qv() |>
    mutate(filt = 'Pass') # |> assign_filters()
describe(summaries$sniffles)

summaries$sniffles2 <- read.csv(
    sprintf('%s/hifi-phasing/summaries/merged.eval.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs(THRESHOLDS) |> bound_qv() |>
    mutate(filt = 'Pass') # |> assign_filters()
describe(summaries$sniffles2)

summaries$illumina_cram <- read.csv(
    sprintf('%s/summaries/327/eval/loo_illumina_cram_v0.17.3.csv.gz', proj_dir),
        sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    group_qvs(THRESHOLDS) |> bound_qv() |>
    assign_filters()
describe(summaries$illumina_cram)

summaries$avail <- mutate(summaries$illumina_loo,
    edit = avail_edit,
    size = avail_size,
    div = edit / size,
    qv = avail_qv, filt = 'Pass') |>
    group_qvs(THRESHOLDS) |>
    bound_qv()
describe(summaries$avail)

local({
    samples <- read_lines('~/Data/data/samples/HPRC.txt')
    summaries$trios <<- read.csv(file.path(eval_dir, 'trio_conc.v0.17.3.csv.gz'),
            sep = '\t', comment = '#') |>
        filter(locus %in% cmrg_loci) |>
        filter(!(indiv %in% samples | father %in% samples | mother %in% samples)) |>
        group_qvs(THRESHOLDS, haps = T) |>
        bound_qv() |>
        mutate(filt = 'Pass')
    cat(sprintf('%s trios loaded\n', length(unique(summaries$trios$indiv))))
})

##########################

plot_order <- matrix(c(
    'illumina_loo', 'Illumina', '',
    'sim_loo', 'Simulated', '',
    'hifi_loo', 'HiFi', '',
    'ont_loo', 'ONT', '',
    'avail', 'Any', 'Avail.',
    'nygc', 'Illumina', '1KGP',
    'sniffles', 'HiFi', 'Sniffles',
    'sniffles2', 'HiFi', 'Sniffles\n+Deepv.',
    'trios', 'Illumina', 'Trio\nconc.',
    'illumina', 'Illumina', '',
    'sim', 'Simulated', '',
    'hifi', 'HiFi', '',
    'ont', 'ONT', ''
), nrow = 3)

aggrs <- lapply(names(summaries),
    function(name) aggr_summaries(summaries[[name]], name)) %>%
    do.call(rbind, .) |>
    filter(key != 'illumina_cram') |>
    mutate(
        facet = case_when(
            grepl('nygc|sniffles', key) ~ 'Other',
            key == 'avail' ~ 'Avail',
            grepl('^(illumina|hifi|ont|sim)_loo$', key) ~ 'Locityper (LOO)',
            key == 'trios' ~ 'Trios',
            grepl('^(illumina|hifi|ont|sim)$', key) ~ 'Locityper (full)',
            ) |>
            factor(levels = c('Locityper (LOO)', 'Avail', 'Other', 'Trios', 'Locityper (full)')),
        x = match(key, plot_order[1,])
    )

(top_half <- ggplot(aggrs) +
    geom_bar(aes(x, frac, fill = qv_cat),
        position = position_stack(reverse = F),
        stat = 'identity', width = 1, color = '#000000C0', linewidth = 0.1) +
    facet_grid2(. ~ facet, space = 'free', scales = 'free', axes = 'y') +
    scale_fill_manual('QV',
        breaks = levels(summaries$illumina$qv_cat),
        values = c(COLORS, 'gray80'),
        drop = F) +
    geom_hline(yintercept = seq(0.2, 0.8, 0.2), linetype = '22',
        color = 'black', alpha = 0.7, linewidth = 0.2) +
    geom_hline(yintercept = c(0, 1),
        color = 'black', alpha = 0.7, linewidth = 0.2) +
    scale_x_continuous('Technology',
        expand = c(0, 0),
        breaks = 1:ncol(plot_order),
        labels = plot_order[2,],
        sec.axis = dup_axis(name = NULL, labels = plot_order[3,]),
        ) +
    scale_y_continuous('Fraction of haplotypes',
        breaks = seq(0, 1, 0.2), expand = c(0, 0),
        ) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        plot.subtitle = element_text(size = 9, margin = margin(), hjust = 0.5),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        # legend.position = 'right',
        # legend.margin = margin(l = -5, r = 0),
        # legend.title = element_text(size = 10, family = 'Overpass SemiBold'),
        # legend.text = element_text(size = 8, family = 'Overpass',
        #     margin = margin(l = 2), hjust = 0.5),
        # legend.key.size = unit(0.8, 'lines'),
        # legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(0.1, 'lines'),
        axis.title.x = element_text(size = 10, margin = margin(t = -4)),
        axis.text.x.bottom = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),
        axis.text.x.top = element_text(size = 8, angle = 90, hjust = 0, vjust = 0.5,
            lineheight = 0.7),
        axis.ticks.x.top = element_blank(),
        axis.ticks.length.x.top = unit(0, 'lines'),
        axis.ticks.length.y = unit(0.15, 'lines'),
        axis.title.y = element_text(size = 11, margin = margin(r = 1), hjust = 0.5),
        axis.text.y = element_text(size = 8),
        plot.background = element_rect(fill = NA, color = NA),
        plot.margin = margin(l = 2, r = 3, t = 2, b = 2)
    ))
ggsave(file.path(plots_dir, 'fig2_top.png'), top_half, bg = 'white',
    width = 6.4, height = 5, dpi = 600, scale = 0.7)

#####################

lost_qv <- rbind(add_column(summaries$hifi_loo, tech = 'HiFi'),
        add_column(summaries$illumina_loo, tech = 'Illumina')) |>
    mutate(tech = factor(tech, levels = c('Illumina', 'HiFi'))) |>
    filter(filt == 'Pass') |>
    bound_qv() |>
    bound_qv(prefix = 'avail_') |>
    mutate(
        round_avail_qv = round(avail_qv),
        corrected_qv = qv + round_avail_qv - avail_qv
    )

MIN_QV <- 17
MAX_QV <- 43

(lost_qv_panel <- ggplot(lost_qv) +
    geom_abline(intercept = c(0, -5, -10), linewidth = 0.3, color = 'firebrick4') +
    geom_violin(aes(round_avail_qv, corrected_qv, group = round_avail_qv), width = 4,
        data = filter(lost_qv, between(round_avail_qv, MIN_QV, MAX_QV)),
        kernel = 'g', adjust = 0.2, bw = 1,
        fill = rgb(59, 61, 126), color = rgb(59, 61, 126),
        linewidth = 0.05, alpha = 0.8,
        ) +
    # geom_point(aes(x = x, y = x), data = data.frame(x = MIN_QV:MAX_QV), color = 'firebrick4') +
    ggh4x::facet_grid2(~ tech, axes = 'y') +
    # geom_point(aes(avail_qv, qv), size = 2, alpha = 0.5, shape = '.') +
    geom_text(aes(x = MIN_QV - 1, y = MAX_QV, label = tech),
        data = data.frame(tech = unique(lost_qv$tech)),
        size = 7.5, family = 'Overpass', hjust = 0, vjust = 1) +
    geom_text(aes(x = MIN_QV - 1.5, y = MIN_QV + loss - 3,
        label = sprintf('%d', loss)),
        data = data.frame(loss = c(0, -5, -10)),
        size = 5, family = 'Overpass', hjust = 0, angle = 18, color = 'firebrick4') +
    scale_x_continuous('Haplotype availabilty (QV)', expand = c(0, 0),
        breaks = seq(20, 40, 10)) +
    scale_y_continuous('Haplotyping accuracy (QV)', expand = c(0, 0)) +    
    coord_cartesian(xlim = c(MIN_QV - 2, MAX_QV + 2), ylim = c(-3, MAX_QV + 1)) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0.7, 'lines')
    )
)
ggsave(file.path(plots_dir, 'fig2_btm.png'), lost_qv_panel, bg = 'white',
    width = 7.8, height = 2.5, dpi = 600, scale = 2)

#################

#################

barplots <- list()
for (key in names(summaries)) {
    barplots[[key]] <- draw_barplot(summaries[[key]])
}

plot_grid(
    barplots$sim_loo + labs(subtitle = 'Simulated data'),
    barplots$hifi_loo + labs(subtitle = 'PacBio HiFi'),
    nrow = 2,
    labels = c('a', 'b'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )
ggsave(sprintf('%s/barplots/supp1ab.png', plots_dir),
    width = 10, height = 13, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/supp1c.png', plots_dir),
    barplots$ont_loo,
    width = 10, height = 6.5, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/avail.png', plots_dir),
    barplots$avail,
    width = 10, height = 6.5, dpi = 600, scale = 0.75)



# ggsave(sprintf('%s/barplots/illumina.png', plots_dir),
#     barplots$illumina,
#     width = 10, height = 7, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/illumina_loo.png', plots_dir),
    barplots$illumina_loo,
    width = 10, height = 7, dpi = 600, scale = 0.75)
# ggsave(sprintf('%s/barplots/avail.png', plots_dir),
#     barplots$avail,
#     width = 10, height = 7, dpi = 600, scale = 0.75)

ggsave(sprintf('%s/barplots/avail.png', plots_dir),
    draw_barplot(summaries$avail, fill_label = 'Availability (QV)'),
    width = 10, height = 6.5, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/nygc.png', plots_dir),
    barplots$nygc,
    width = 10, height = 6.5, dpi = 600, scale = 0.75)

plot_grid(
    barplots$avail + labs(subtitle = 'Locityper, Haplotype availability'),
    barplots$sim_loo + labs(subtitle = 'Sim'),
    nrow = 2,
    labels = c('a', 'b'),
    label_fontfamily = 'Overpass',
    label_y = 1.015
    )

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
# draw_barplot(counts$trios, fill_label = 'Concordance (QV)')

##################

filter(summaries$illumina_loo, filt == 'Pass') |>
    group_by(locus) |>
    filter(mean(qv_num <= 2) >= 0.9) |>
    with(length(unique(locus)))

describe_several(summaries, 'loo')
describe(summaries$illumina_cram)
describe(summaries$avail)
describe

describe_several(summaries, '_loo|avail')
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

#####

.summary_to_row <- function(summary, key, overall_total = NA) {
    tab <- count(summary, qv_cat, name = 'count') |>
        mutate(perc = 100 * count / sum(count))
    row <- pivot_wider(tab, names_from = 'qv_cat', values_from = c('count', 'perc')) |>
        mutate(key = key, mean = mean(summary$qv), median = median(summary$qv),
            total = nrow(summary), total_perc = 100 * total / overall_total)
    cols <- c('key', 'mean', 'median', 'total', 'total_perc')
    for (cat in tab$qv_cat) {
        cols <- c(cols, sprintf('%s_%s', c('count', 'perc'), cat))
    }
    row[cols]
}
summary_to_row <- function(summary, key) {
    if (any(summary$qv_cat == 'No call')) {
        rbind(
            mutate(summary, qv_cat = qv_cat0) |> .summary_to_row(key),
            filter(summary, qv_cat != 'No call') |> .summary_to_row(key, nrow(summary))
        )
    } else {
        .summary_to_row(summary, key)
    }
}
tab <- lapply(names(summaries),
    function(name) summary_to_row(summaries[[name]], name)) %>%
    do.call(rbind, .) |>
    mutate(key = factor(key, levels = c(
        'illumina_loo', 'sim_loo', 'hifi_loo', 'ont_loo', 'illumina_cram',
        'avail',
        'nygc', 'sniffles', 'sniffles2',
        'trios',
        'illumina', 'sim', 'hifi', 'ont'
        ))) |>
    arrange(key, is.na(total_perc))
write.table(tab, file.path(plots_dir, 'barplots.csv'), row.names = F, quote = F, sep = '\t')
