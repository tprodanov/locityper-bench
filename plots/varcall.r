if (!exists('.IMPORTED_COMMON')) { source('common.r') }
library(ggh4x)
library(see)

qual_threshs <- data.frame(
    tool = c(rep(c('locityper', 'locityper_loo', 'pangenie'), each = 2), 'nygc'),
    threshold = c(0, 1, 0, 1, 0, 10, 0),
    filt = c(rep(c('Any', 'High'), 3), 'High')
)
acc <- read.csv(file.path(proj_dir, 'comparison/VCF/combined_evals.csv.gz'),
    sep = '\t', comment = '#') |>
    filter(locus %in% cmrg_loci) |>
    mutate(fn = replace_na(fn, 0)) |>
    left_join(qual_threshs, by = c('tool', 'threshold')) |>
    filter(!is.na(filt))
sum_acc <- group_by(acc, locus, tool, type, filt) |>
    summarize(
        tp_base = sum(tp_base), tp_call = sum(tp_call),
        fp = sum(fp), fn = sum(fn), .groups = 'keep') |>
    ungroup() |>
    mutate(
        precision = ifelse(tp_call + fp > 0, tp_call / (tp_call + fp), 1),
        recall = ifelse(tp_base + fn > 0, tp_base / (tp_base + fn), 1),
        f1 = ifelse(precision + recall > 0,
            2 * precision * recall / (precision + recall), 0))
wsum_acc <- select(sum_acc, locus, tool, type, filt, precision, recall, f1) |>
    pivot_longer(c(precision, recall, f1), names_to = 'metric') |>
    filter(tool != 'locityper_loo') |>
    mutate(
        type = recode_factor(type,
            'snps' = 'SNPs',
            'indels' = 'Indels  &  SVs',
            'all' = 'All'),
        metric = recode_factor(metric,
            'precision' = 'Precision',
            'recall' = 'Recall',
            'f1' = 'F₁  score'),
        filt = factor(filt, levels = c('Any', 'High')),
        tool = recode_factor(tool,
            'locityper' = 'Locityper',
            'pangenie' = 'Pangenie',
            'nygc' = '1KGP'))

colors <- c(
    'Locityper' = rgb(255, 221, 35),
    'Pangenie' = rgb(104, 193, 189),
    '1KGP' = rgb(121, 109, 180)
    )

ggplot(filter(wsum_acc, filt == 'High'),
        aes(tool, value, fill = tool)) +
    geom_violin(
        linewidth = 0.3,
        width = 1.05,
        kernel = 'g',
        bw = 0.3,
        adjust = 0.3,
        draw_quantiles = 0.5,
        alpha = 0.9) +
    facet_nested(vars(type), vars(metric)) +
    scale_x_discrete(NULL) +
    scale_y_continuous('Metric value', breaks = c(0, 1), expand = expansion(add = 0.05)) +
    scale_fill_manual('Call set ', values = colors) +
    scale_alpha_manual('    Quality ', values = c(0.3, 1), breaks = c('Any', 'High')) +
    guides(
        fill = guide_legend(order = 1),
        alpha = guide_legend(order = 2, override.aes = list(fill = 'gray30'))) +
    guides(
        fill = guide_legend(order = 1),
        alpha = guide_legend(order = 2, override.aes = list(fill = 'gray30'))) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        strip.background = element_rect(color = NA, fill = 'gray90'),
        strip.text = element_text(family = 'Overpass Semibold', margin = margin(2, 2, 2, 2)),
        strip.placement = 'outside',
        panel.border = element_rect(color = NA),
        panel.spacing.x = unit(0.8, 'lines'),
        panel.spacing.y = unit(0.2, 'lines'),

        plot.margin = margin(2, 2, 2, 2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = 'gray80', linewidth = 0.4),
        panel.grid.minor.y = element_line(color = 'gray80', linewidth = 0.3),
        legend.position = 'none',
        axis.text.x = element_text(size = 8),
    )
ggsave(file.path(plots_dir, 'varcall.png'),
    width = 8, height = 4.5, dpi = 600, scale = 0.8)

filter(wsum_acc, filt == 'High') |>
    group_by(tool, metric, type) |>
    summarize(mean = round(mean(value), 3), median = round(median(value), 3)) |>
    ungroup() |>
    View()
