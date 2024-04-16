if (!exists('.IMPORTED_COMMON')) { source('common.r') }

MAIN_FIG <- F

extend_accuracy <- function(accuracies) {
    if (MAIN_FIG) {
        max_acc <- 4
        max_acc_str <- 'Full match'
        mutate(accuracies,
            accuracy2 = case_when(
                startsWith(base, '<') & startsWith(pred, '<') ~ max_acc_str,
                accuracy >= max_acc ~ max_acc_str,
                accuracy > 0 & !is.na(base_length) & accuracy == base_length ~ max_acc_str,
                T ~ as.character(accuracy)
            ) |> factor(levels = c(as.numeric((0:max_acc) - 1), max_acc_str)))
    } else {
        mutate(accuracies,
            accuracy2 = case_when(
                startsWith(base, '<') & startsWith(pred, '<') ~ 'Deletion found',
                startsWith(base, '<') ~ 'Extra copy',
                startsWith(pred, '<') ~ 'Missing copy',
                accuracy >= 2 ~ 'Correct protein',
                !is.na(avail) & !avail ~ 'Unavailable protein',
                accuracy < 2 ~ 'Incorrect protein'
            ) |> factor(levels = c('Extra copy', 'Missing copy',
                'Incorrect protein', 'Unavailable protein',
                'Correct protein', 'Deletion found')))
    }
}

traffic <- ggthemes::tableau_color_pal('Traffic')(9)
hue_circle <- ggthemes::tableau_color_pal('Hue Circle')(19)
superfishel <- ggthemes::tableau_color_pal('Superfishel Stone')(10)

if (MAIN_FIG) {
    colors3 <- c(
        '0' = traffic[7],
        '1' = hue_circle[9],
        '2' = hue_circle[7],
        '3' = hue_circle[6],
        'Full match' = '#437b29'
    )
    fill_name <- 'Correctly predicted fields'
    lw <- 0.1
} else {
    colors3 <- c(
        'Extra copy' = '#606060', # superfishel[10],
        'Missing copy' = superfishel[2],
        'Incorrect protein' = superfishel[3],
        'Unavailable protein' = superfishel[8],
        'Correct protein' = superfishel[4],
        'Deletion found' = superfishel[1]
    )
    fill_name <- NULL
    lw <- 0.2
}

hla_get_counts <- function(acc) {
    totals <- count(acc, tool, gene, name = 'total')
    count(acc, tool, gene, accuracy2) |>
        arrange(tool, desc(accuracy2)) |>
        left_join(totals, by = c('tool', 'gene')) |>
        mutate(frac = n / total) |>
        group_by(tool, gene) |>
        mutate(cum_frac = cumsum(frac)) |>
        ungroup()
}

draw_accuracy <- function(counts, nrows, strip.face = 'italic') {
    ggplot(counts) +
        geom_bar(aes(tool, frac, fill = accuracy2), stat = 'identity',
            color = 'black', linewidth = 0.1, width = 1) +
        geom_hline(yintercept = seq(0.25, 0.75, 0.25),
            linetype = '14', linewidth = 0.3, color = 'gray10') +
        facet_wrap(~ gene, nrow = nrows) +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous(
            ifelse(nrows == 1, 'Frac. of haplotypes     ', 'Fraction of haplotypes'),
            breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
        scale_fill_manual(fill_name, values = colors3) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            legend.title = element_text(margin = margin(r = 5)),
            legend.position = 'top',
            legend.margin = margin(b = -5),
            legend.key.size = unit(0.8, 'lines'),
            legend.text = element_text(margin = margin(l = -1, r = 4)),
            axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
            strip.text = element_text(face = strip.face,
                size = 9, margin = margin(t = 1, b = 1)),
            strip.background = element_rect(fill = 'gray95', color = NA),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = NA, color = NA),
            panel.border = element_rect(fill = NA, color = NA),
            panel.spacing.y = unit(0.2, 'lines'),
        )
}

#############

hla_dir <- file.path(proj_dir, 'comparison/HLA')
loci <- read.csv(file.path(hla_dir, 'loci.txt'), sep = '\t', header = F)
samples <- readLines(file.path(proj_dir, 'samples/illumina_40.txt'))

acc <- rbind(
    read.csv(file.path(hla_dir, 'eval/locityper.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'Locityper ●'),
    read.csv(file.path(hla_dir, 'eval/locityper_loo.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'Locityper ○'),
    read.csv(file.path(hla_dir, 'eval/t1k.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'T1K')) |>
    filter(gene %in% loci$V1 & sample %in% samples) |>
    mutate(
        tool = factor(tool, levels = unique(tool)),
        ) |>
    extend_accuracy()

counts1 <- filter(acc, !grepl('^KIR', gene)) |>
    mutate(gene = sub('^HLA-', '-', gene)) |>
    hla_get_counts()
counts2 <- filter(acc, grepl('^KIR', gene)) |>
    mutate(gene = sub('^KIR', '', gene)) |>
    hla_get_counts()
counts3 <- mutate(acc,
        gene = factor(ifelse(grepl('^KIR', gene), 'KIR', 'MHC'), levels = c('MHC', 'KIR'))) |>
    hla_get_counts()
plot_grid(
    draw_accuracy(counts1, 3),
    plot_grid(
        draw_accuracy(counts2, 1) + theme(legend.position = 'none'),
        draw_accuracy(counts3, 1, 'bold') + theme(legend.position = 'none'),
        rel_widths = c(0.75, 0.25),
        labels = c('b', 'c'),
        label_fontfamily = 'Overpass',
        label_y = 1.06,
        label_x = c(-0.007, -0.017)
        ),
    rel_heights = c(0.67, 0.33),
    nrow = 2,
    labels = 'a',
    label_x = -0.003,
    label_y = 1.01,
    label_fontfamily = 'Overpass'
    )
ggsave(sprintf('%s/HLA_%s.png', plots_dir, ifelse(MAIN_FIG, 'main', 'supp')),
    width = 10, height = 7, dpi = 600, scale = 0.8)

#############

filter(counts3, gene == 'MHC')
filter(counts3, gene == 'KIR')

filter(counts3, accuracy2 != 'Correct group' & accuracy2 != 'Deletion found') |>
    group_by(tool, gene) |>
    reframe(accuracy2 = accuracy2, m = round(n / sum(n), 3)) |>
    ungroup() |>
    as.data.frame()
