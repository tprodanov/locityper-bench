if (!exists('.IMPORTED_COMMON')) { source('common.r') }

time <- read.csv('~/Data/proj/locityper/mutated/time.csv', sep = '\t') |>
    mutate(step = ifelse(step == 'init', 'recruit', step),
        step = recode_factor(step,
            'recruit' = 'Read recruitment', 'map' = 'Read mapping', 'genotype' = 'Genotyping',
            'total' = 'Total'),
        minutes = seconds / 60)
sum_time <- aggregate(minutes ~ size + sample + step, time, sum) |>
    mutate(group = paste0(sample, step))
aver_time <- aggregate(minutes ~ size + step, sum_time, mean)
norm_fct <- filter(aver_time, step == 'Total' & size == 1)$minutes

(panel_a <- ggplot(sum_time, aes(size, minutes, color = step)) +
    # geom_point(alpha = 0.5) +
    geom_line(aes(group = group), alpha = 0.5) +
    geom_point(data = aver_time, size = 2) +
    geom_line(data = aver_time, linewidth = 1.5) +
    scale_x_continuous('Relative reference panel size',
        breaks = 1:10, labels = function(x) sprintf('×%s', x), expand = expansion(mult = 0.01)) +
    scale_y_continuous('Runtime (minutes)', limits = c(0, NA), expand = expansion(mult = 0),
        sec.axis = sec_axis(~ . / norm_fct, name = 'Normalized runtime',
            breaks = 1:10,
            labels = function(x) sprintf('×%s', x))
        ) +
    scale_color_manual(NULL,
        values = c('#006e00', '#b80058', '#008cf9', '#ebac23') # https://tsitsul.in/blog/coloropt/
        ) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        legend.position = 'top',
        legend.margin = margin(b = -10),
        legend.text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.y.right = unit(0, 'lines')
    ))

mem <- read.csv('~/Data/proj/locityper/mutated/mem.csv', sep = '\t') |>
    mutate(gb = memory / 1024)
aver_mem <- aggregate(gb ~ size, mem, mean)
norm_fct2 <- filter(aver_mem, size == 1)$gb

(panel_b <- ggplot(mem, aes(size, gb)) +
    geom_line(aes(group = sample), alpha = 0.5) +
    geom_point(data = aver_mem, size = 2) +
    geom_line(data = aver_mem, linewidth = 1.5) +
    scale_x_continuous('Relative reference panel size',
        breaks = 1:10, labels = function(x) sprintf('×%s', x), expand = expansion(mult = 0.01)) +
    scale_y_continuous('Peak memory (Gb)', limits = c(0, norm_fct2 * 5),
        expand = expansion(mult = 0),
        sec.axis = sec_axis(~ . / norm_fct2, name = 'Normalized peak memory',
            breaks = 1:10,
            labels = function(x) sprintf('×%s', x))
        ) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        legend.position = 'top',
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length.y.right = unit(0, 'lines')
    ))
cowplot::plot_grid(panel_a, panel_b,
    nrow = 2,
    labels = letters,
    label_fontfamily = 'Overpass',
    rel_heights = c(0.6, 0.4),
    align = 'v')
ggsave(file.path(plots_dir, 'time.png'), width = 8, height = 10, dpi = 600, scale = 0.7)
