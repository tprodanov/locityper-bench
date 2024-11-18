if (!exists('.IMPORTED_COMMON')) { source('common.r') }

zissou <- wesanderson::wes_palette("Zissou1", 20, type = "continuous")
zissou100 <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
scales::show_col(zissou100, F)
THRESHOLDS <- c(0, 17, 23, 33, 43)
COLORS <- zissou100[c(1, 11, 35, 65, 92)]
# scales::show_col(COLORS, F)
stopifnot(length(THRESHOLDS) == length(COLORS))

draw_barplot <- function(counts, panel_width = 71, fill_label = 'Accuracy (QV)')
{
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
describe <- function(counts) {
    cats <- levels(counts$qv_cat)
    n <- sum(counts$count)
    cat(sprintf('%-12s  %-20s  %-20s  %-20s\n',
        '', 'this accuracy', 'this or more acc.', 'this or less acc.'))
    for (i in 1:length(cats)) {
        cat(sprintf('    %-8s  %20s  %20s  %20s\n',
            ascii_qvs(cats[i]),
            percentage(with(counts, sum((qv_num == i) * count)), n),
            percentage(with(counts, sum((qv_num <= i) * count)), n),
            percentage(with(counts, sum((qv_num >= i) * count)), n)
            ))
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

vers <- c('v0.17.0') #, 'v0.14.2-4')
summaries <- list()
counts <- list()
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
            counts[[key]] <- get_counts(summaries[[key]])
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
counts$nygc <- get_counts(summaries$nygc)
describe(counts$nygc)

summaries$avail <- mutate(summaries$illumina_loo,
    edit = avail_edit,
    size = avail_size,
    div = edit / size,
    qv = avail_qv, filt = 'Pass') |>
    group_qvs(THRESHOLDS) |>
    bound_qv()
counts$avail <- get_counts(summaries$avail)
describe(counts$avail)

local({
    summaries$trios <<- read.csv(file.path(eval_dir, 'trio_conc.csv.gz'),
            sep = '\t', comment = '#') |>
        filter(locus %in% cmrg_loci) |>
        group_qvs(THRESHOLDS, haps = T) |>
        bound_qv() |>
        mutate(filt = 'Pass')
    counts$trios <<- get_counts(summaries$trios)
})

barplots <- list()
for (key in names(counts)) {
    barplots[[key]] <- draw_barplot(counts[[key]])
}

##########################

aggrs <- lapply(names(summaries),
    function(name) aggr_summaries(summaries[[name]], name)) %>%
    do.call(rbind, .) |>
    mutate(
        panel = case_when(
            key == 'avail' ~ 'Avail.',
            key == 'trios' ~ 'Trios',
            T ~ 'Genotyping'),
        facet = case_when(
            grepl('^(illumina|hifi|ont|sim)_loo$', key) ~ 'Locityper (LOO)',
            grepl('^(illumina|hifi|ont|sim)$', key) ~ 'Locityper (full)',
            key == 'nygc' ~ '1KGP',
            T ~ ' ') |>
            factor(levels = c('Locityper (LOO)', '1KGP', 'Locityper (full)', ' ')),
        label = case_when(
            grepl('illumina', key) ~ 'Illumina',
            grepl('sim', key) ~ 'Simulated',
            grepl('hifi', key) ~ 'HiFi',
            grepl('ont', key) ~ 'ONT',
            key == 'avail' ~ 'Any',
            key == 'nygc' ~ 'Illumina',
            key == 'trios' ~ 'Illumina') |>
            factor(levels = c('Illumina', 'Simulated', 'HiFi', 'ONT', 'Any'))
    )

aggr_plots <- list()
for (title in unique(aggrs$panel)) {
    aggr_plots[[title]] <- filter(aggrs, panel == title) |>
        ggplot() +
        geom_bar(aes(label, frac, fill = qv_cat),
            position = position_stack(reverse = F),
            stat = 'identity', width = 1, color = '#000000C0', linewidth = 0.1) +
        # {
        #     if (title == 'Genotyping') { facet_grid(. ~ facet, space = 'free', scales = 'free') }
        #     else { NULL }
        # } +
        facet_grid(. ~ facet, space = 'free', scales = 'free') +
        scale_fill_manual('QV',
            breaks = c(levels(summaries$illumina$qv_cat), 'No call'),
            values = c(COLORS, 'gray80'),
            drop = F) +
        geom_hline(yintercept = seq(0.2, 0.8, 0.2), linetype = '22',
            color = 'black', alpha = 0.7, linewidth = 0.2) +
        geom_hline(yintercept = c(0, 1),
            color = 'black', alpha = 0.7, linewidth = 0.2) +
        scale_x_discrete(if (title == 'Genotyping') { 'Technology' } else { NULL },
            expand = c(0, 0)) +
        scale_y_continuous('Fraction of haplotypes',
            breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
        ggtitle(NULL, subtitle = title) +
        # guides(fill = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            plot.subtitle = element_text(size = 9, margin = margin(), hjust = 0.5),
            panel.background = element_rect(fill = NA),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'right',
            legend.margin = margin(l = 10, r = 5),
            legend.title = element_text(size = 10, family = 'Overpass SemiBold'),
            legend.text = element_text(size = 8, family = 'Overpass',
                margin = margin(l = 2), hjust = 0.5),
            legend.key.size = unit(0.8, 'lines'),
            legend.background = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 8, margin = margin(t = 2, b = 2), hjust = 0.5),
            panel.spacing.x = unit(0.3, 'lines'),
            axis.title.x = element_text(size = 10, margin = margin()),
            axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
            plot.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(l = 3, r = 3, t = 2, b = 2)
        )
}
aggr_plots$Genotyping

# (top_half <- plot_grid(
#     plot_grid(
#         aggr_plots$`Avail.` + theme(legend.position = 'none') + labs(x = 'NULL'),
#         aggr_plots$Genotyping + theme(legend.position = 'none', axis.title.y = element_blank()),
#         aggr_plots$Trios + theme(legend.position = 'none', axis.title.y = element_blank()),
#         nrow = 1,
#         rel_widths = c(0.13, 0.45, 0.1),
#         align = 'h'
#     ),
#     get_plot_component(aggr_plots$Genotyping, 'guide-box', return_all = T)[[1]],
#     nrow = 1,
#     rel_widths = c(0.87, 0.13)
#     ))
(top_half <- plot_grid(
    aggr_plots$`Avail.` + theme(legend.position = 'none') + labs(x = 'NULL'),
    aggr_plots$Genotyping + theme(legend.position = 'none', axis.title.y = element_blank()),
    aggr_plots$Trios + theme(legend.position = 'none', axis.title.y = element_blank()),
    nrow = 1,
    rel_widths = c(0.13, 0.45, 0.1),
    align = 'h'
    ))
ggsave(file.path(plots_dir, 'fig2_top.png'), top_half, bg = 'white',
    width = 6, height = 5, dpi = 600, scale = 0.7)

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
        #strip.text = element_text(size = 20, hjust = 0.5, margin = margin(t = 2, b = -15)),
        strip.text = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        # panel.border = element_rect(color = 'gray80', linewidth = 0.3),
        #axis.line = element_line(color = 'gray30', linewidth = 0.3),
        panel.spacing.x = unit(0.7, 'lines')
    )
)
ggsave(file.path(plots_dir, 'fig2_btm.png'), lost_qv_panel, bg = 'white',
    width = 7.8, height = 2.5, dpi = 600, scale = 2)

###############

# ext_counts <- rbind(
#     mutate(counts$illumina_loo, ty = 'accuracy'),
#     mutate(counts$avail, ty = 'avail')
# )
# n_loci <- length(levels(ext_counts$locus))
# panels <- 5
# panel_size <- ceiling(n_loci / panels)
# 
# width <- 0.9
# half_width <- width / 2
# ext_counts <- mutate(ext_counts,
#     locus_ix = as.numeric(locus),
#     x = locus_ix + ifelse(ty == 'avail', -1, 1) * half_width / 2,
#     panel = (locus_ix - 1) %/% panel_size
# )
# 
# ggplot(ext_counts) +
#     geom_bar(aes(x, frac, fill = qv_cat),
#         stat = 'identity', position = position_stack(reverse = F),
#         width = half_width, color = '#000000C0', linewidth = 0.1) +
#     # geom_hline(yintercept = seq(0, 1, 0.25), linetype = '22',
#     #         color = '#000000C0', linewidth = 0.2) +
#     facet_wrap(~ panel, ncol = 1, scales = 'free_x') +
#     scale_x_continuous(NULL, expand = c(0, 0),
#         breaks = 1:n_loci, labels = levels(ext_counts$locus)) +
#     scale_y_continuous('Fraction of haplotypes', expand = c(0, 0)) +
#     scale_fill_manual('Accuracy/Availability (QV)',
#         breaks = levels(ext_counts$qv_cat),
#         values = c(COLORS, 'gray80')) +
#     guides(fill = guide_legend(nrow = 1)) +
#     theme_bw() +
#     theme(
#         text = element_text(family = 'Overpass'),
#         panel.grid = element_blank(),
#         panel.background = element_rect(fill = 'white'),
#         panel.border = element_blank(),
#         legend.position = 'bottom',
#         legend.margin = margin(t = -10),
#         legend.title = element_text(size = 11, family = 'Overpass SemiBold'),
#         legend.text = element_text(family = 'Overpass',
#             margin = margin(l = 2, r = 3)),
#         legend.key.size = unit(0.9, 'lines'),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         panel.spacing.y = unit(0.1, 'lines'),
#         axis.text.y = element_text(size = 8),
#         axis.text.x = element_text(
#             size = 6.2, face = 'italic', angle = 45, vjust = 1, hjust = 1)
#     )
# ggsave(sprintf('%s/barplots/loo+avail.png', plots_dir),
#     width = 10, height = 9, dpi = 600, scale = 0.75)
# 
# loo_ext <- summaries$illumina_loo
# n_loci <- length(levels(loo_ext$locus))
# panels <- 4
# panel_size <- ceiling(n_loci / panels)
# 
# loo_ext <- mutate(loo_ext,
#     locus_ix = as.numeric(locus),
#     panel = (locus_ix - 1) %/% panel_size
# )
# loo_ext <- mutate(summaries$illumina_loo

###############

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
    draw_barplot(counts$avail, fill_label = 'Availability (QV)'),
    width = 10, height = 6.5, dpi = 600, scale = 0.75)
ggsave(sprintf('%s/barplots/nygc.png', plots_dir),
    draw_barplot(counts$nygc),
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
