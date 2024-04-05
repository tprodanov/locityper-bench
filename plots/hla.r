if (!exists('.IMPORTED_COMMON')) { source('common.r') }

SHOW_DELS <- F

extend_accuracy <- function(accuracies) {
    if (SHOW_DELS) {
        mutate(accuracies,
            accuracy2 = case_when(
                startsWith(base, '<') & startsWith(pred, '<') ~ 'Deletion found',
                startsWith(base, '<') ~ 'Extra copy',
                startsWith(pred, '<') ~ 'Missing copy',
                accuracy > 0 ~ 'Correct group',
                !is.na(avail) & !avail ~ 'Unavailable group',
                accuracy == 0 ~ 'Incorrect group'
            ) |> factor(levels = c('Extra copy', 'Missing copy',
                'Incorrect group', 'Unavailable group', 'Correct group', 'Deletion found')))
    } else {
        mutate(accuracies,
            accuracy2 = case_when(
                startsWith(base, '<') & startsWith(pred, '<') ~ 'Deletion found',
                T ~ as.character(accuracy)
            ) |> factor(levels = c('0', '1', '2', '3', 'Deletion found')))
    }
}

traffic <- ggthemes::tableau_color_pal('Traffic')(9)
hue_circle <- ggthemes::tableau_color_pal('Hue Circle')(19)
superfishel <- ggthemes::tableau_color_pal('Superfishel Stone')(10)

if (SHOW_DELS) {
    colors3 <- c(
        'Extra copy' = superfishel[10],
        'Missing copy' = superfishel[2],
        'Incorrect group' = superfishel[3],
        'Unavailable group' = superfishel[8],
        'Correct group' = superfishel[4],
        'Deletion found' = superfishel[1]
    )
    fill_name <- NULL
    lw <- 0.2
} else {
    colors3 <- c(
        'Missed copy' = traffic[1],
        'Missed DEL' = hue_circle[17],
        '0' = traffic[7],
        '1' = hue_circle[9],
        '2' = hue_circle[7],
        '3' = hue_circle[6],
        'Deletion found' = rgb(0, 139, 176)
    )
    fill_name <- 'Accuracy level'
    lw <- 0.1
}


draw_accuracy <- function(acc, nrows, strip.face = 'italic') {
    totals <- count(acc, tool, gene, name = 'total')
    counts <- count(acc, tool, gene, accuracy2) |>
        left_join(totals, by = c('tool', 'gene')) |>
        mutate(frac = n / total)

    ggplot(counts) +
        geom_bar(aes(tool, frac, fill = accuracy2), stat = 'identity',
            color = 'black', linewidth = 0.1, width = 1) +
        geom_hline(yintercept = seq(0.25, 0.75, 0.25),
            linetype = '14', linewidth = 0.3, color = 'gray10') +
        facet_wrap(~ gene, nrow = nrows) +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous(
            ifelse(nrows == 1, 'Frac. of haplotypes', 'Fraction of haplotypes'),
            breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
        scale_fill_manual(fill_name, values = colors3) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            legend.position = 'top',
            legend.margin = margin(b = -5),
            legend.key.size = unit(0.85, 'lines'),
            axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1),
            strip.text = element_text(face = strip.face,
                size = 9, margin = margin(t = 1, b = 1)),
            strip.background = element_rect(fill = 'gray95', color = NA),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = NA, color = NA),
            panel.border = element_rect(fill = NA, color = NA),
        )
}

#############

hla_dir <- file.path(proj_dir, 'comparison/HLA')
loci <- read.csv(file.path(hla_dir, 'loci.txt'), sep = '\t', header = F)
samples <- readLines(file.path(proj_dir, 'samples/illumina_40.txt'))

acc <- rbind(
    read.csv(file.path(hla_dir, 'eval/locityper.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'Locityper'),
    read.csv(file.path(hla_dir, 'eval/locityper_loo.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'Locityper✳'),
    read.csv(file.path(hla_dir, 'eval/t1k.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'T1K')) |>
    filter(gene %in% loci$V1 & sample %in% samples) |>
    mutate(
        accuracy = pmin(accuracy, 3),
        tool = factor(tool, levels = unique(tool)),
        ) |>
    extend_accuracy()


g1 <- filter(acc, !grepl('^KIR', gene)) |>
    mutate(gene = sub('^HLA-', '-', gene)) |>
    draw_accuracy(3)
g2 <- (filter(acc, grepl('^KIR', gene)) |>
    mutate(gene = sub('^KIR', '', gene)) |>
    draw_accuracy(1)) +
    theme(legend.position = 'none')
g3 <- (mutate(acc,
        gene = factor(ifelse(grepl('^KIR', gene), 'KIR', 'MHC'), levels = c('MHC', 'KIR'))) |>
    draw_accuracy(1, 'bold')) +
    theme(legend.position = 'none')
plot_grid(
    g1,
    plot_grid(
        g2, g3,
        rel_widths = c(0.75, 0.25),
        labels = c('b', 'c'),
        label_fontfamily = 'Overpass',
        label_y = 1.06,
        label_x = c(-0.007, -0.017)
        ),
    rel_heights = c(0.65, 0.35),
    nrow = 2,
    labels = 'a',
    label_x = -0.003,
    label_fontfamily = 'Overpass'
    )

ggsave(sprintf('%s/HLA_%s.png', plots_dir, ifelse(SHOW_DELS, 'type', 'acc')),
    width = 10, height = 7, dpi = 600, scale = 0.85)

