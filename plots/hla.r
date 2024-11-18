if (!exists('.IMPORTED_COMMON')) { source('common.r') }
library(ggtext)

extend_accuracy <- function(accuracies) {
    max_acc <- 4
    max_acc_str <- 'Full match'
    mutate(accuracies,
        group1 = case_when(
            startsWith(base, '<') & startsWith(pred, '<') ~ max_acc_str,
            accuracy >= max_acc ~ max_acc_str,
            accuracy > 0 & !is.na(base_length) & accuracy == base_length ~ max_acc_str,
            T ~ as.character(accuracy)
        ) |> factor(levels = c(as.numeric((0:max_acc) - 1), max_acc_str)),
        group2 = case_when(
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

colors <- local({
    traffic <- ggthemes::tableau_color_pal('Traffic')(9)
    hue_circle <- ggthemes::tableau_color_pal('Hue Circle')(19)
    superfishel <- ggthemes::tableau_color_pal('Superfishel Stone')(10)
    list(
        c(
            '0' = traffic[7],
            '1' = hue_circle[9],
            '2' = hue_circle[7],
            '3' = hue_circle[6],
            'Full match' = '#437b29'
        ),
        c(
            'Extra copy' = '#606060', # superfishel[10],
            'Missing copy' = superfishel[2],
            'Incorrect protein' = superfishel[3],
            'Unavailable protein' = superfishel[8],
            'Correct protein' = superfishel[4],
            'Deletion found' = superfishel[1]
        ))
})
fill_names <- list('Correctly predicted fields', NULL)

hla_get_counts <- function(acc, preset) {
    totals <- count(acc, tool, gene, name = 'total')
    rename(acc, group = sprintf('group%d', preset)) |>
        count(tool, gene, group) |>
        arrange(tool, desc(group)) |>
        left_join(totals, by = c('tool', 'gene')) |>
        mutate(frac = n / total) |>
        group_by(tool, gene) |>
        mutate(cum_frac = cumsum(frac)) |>
        ungroup()
}

draw_accuracy <- function(acc, preset, nrows) {
    counts <- hla_get_counts(acc, preset)
    ggplot(counts) +
        geom_bar(aes(tool, frac, fill = group), stat = 'identity',
            color = 'black', linewidth = 0.15, width = 1) +
        geom_hline(yintercept = seq(0.25, 0.75, 0.25),
            linetype = '14', linewidth = 0.3, color = 'gray10') +
        facet_wrap(~ gene, nrow = nrows) +
        scale_x_discrete(NULL, expand = c(0, 0)) +
        scale_y_continuous('Fraction of haplotypes',
            breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
        scale_fill_manual(fill_names[[preset]], values = colors[[preset]]) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(
            text = element_text(family = 'Overpass'),
            legend.position = 'top',
            legend.margin = margin(b = -5),
            legend.key.size = unit(0.8, 'lines'),
            legend.title = element_text(margin = margin(r = 5)),
            legend.text = element_text(margin = margin(l = 3, r = 5)),
            axis.text.x =  element_markdown(size = 7, angle = 90, hjust = 1, vjust = 0.5),
            strip.text = element_markdown(size = 9, margin = margin(t = 1, b = 1)),
            strip.background = element_rect(fill = 'gray95', color = NA),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = NA, color = NA),
            panel.border = element_rect(fill = NA, color = NA),
            panel.spacing.x = unit(0.2, 'lines'),
            panel.spacing.y = unit(0.2, 'lines'),
        )
}

full_fig <- function(acc1, acc2, preset) {
    plot_grid(
        draw_accuracy(acc1, preset, 3),
        plot_grid(
            NULL,
            draw_accuracy(acc2, preset, 2) + theme(legend.position = 'none'),
            NULL,
            nrow = 1,
            rel_widths = c(0.12, 0.78, 0.10)
        ),
        rel_heights = c(0.6, 0.4),
        nrow = 2,
        labels = letters,
        label_x = c(0.01, 0.123),
        label_y = 1.01,
        label_fontfamily = 'Overpass')
}

#############

hla_dir <- file.path(proj_dir, 'MHC')
loci <- read.csv(file.path(hla_dir, 'loci1.txt'), sep = ' ', header = F)$V1
samples <- readLines(file.path(proj_dir, 'samples/illumina_40.txt'))

acc <- rbind(
    read.csv(file.path(hla_dir, 'eval/locityper_full.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'F'),
    read.csv(file.path(hla_dir, 'eval/locityper_loo.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'L'),
    read.csv(file.path(hla_dir, 'eval/locityper_loow.csv'), sep = '\t', comment = '#') |>
        mutate(tool = 'W'),
    read.csv(file.path(hla_dir, 'eval/t1k.csv'), sep = '\t', comment = '#') |>
        mutate(tool = '**T1K**')) |>
    filter(gene %in% loci & sample %in% samples) |>
    mutate(
        tool = factor(tool, levels = unique(tool)),
        ) |>
    extend_accuracy()

acc1 <- filter(acc, !grepl('^KIR', gene)) |>
    mutate(gene = sub('^HLA-', '-', gene) %>% sprintf('*%s*', .) |> factor()) %>%
    rbind(., mutate(., gene = '**Average**'))
acc2 <- filter(acc, grepl('KIR', gene)) |>
    mutate(gene = sub('^KIR', '', gene) %>% sprintf('*%s*', .) |> factor()) %>%
    rbind(., mutate(., gene = '**Average**'))

(main_fig <- full_fig(acc1, acc2, 1))
ggsave(file.path(plots_dir, 'HLA_main.png'), main_fig, bg = 'white',
    width = 10, height = 8, dpi = 600, scale = 0.8)

(supp_fig <- full_fig(acc1, acc2, 2))
ggsave(file.path(plots_dir, 'HLA_supp.png'), supp_fig, bg = 'white',
    width = 10, height = 8, dpi = 600, scale = 0.8)



#############

acc3 <- mutate(acc, gene0 = gene, gene = ifelse(grepl('^KIR', gene0), 'KIR', 'MHC'))
select(acc3, gene0, gene) |> unique() |> count(gene)

counts_a <- hla_get_counts(acc3, 1)
filter(counts_a, gene == 'MHC')
filter(counts_a, gene == 'KIR')

counts_b <- hla_get_counts(acc3, 2)
select(counts_b, !c(total, cum_frac)) |> arrange(tool, gene, group) |> print(n = 100)
filter(counts_b, group != 'Correct protein' & group != 'Deletion found') |>
    group_by(tool, gene) |>
    reframe(group = group, frac = round(100 * n / sum(n), 1)) |>
    ungroup() |>
    as.data.frame()

group_by(counts_b, tool, gene) |>
    reframe(group = group, perc = round(100 * n / sum(n), 1)) |>
    ungroup() |>
    filter(group != 'Correct protein' & group != 'Deletion found') |>
    as.data.frame()

filter(counts_b, group != 'Missing copy' & group != 'Unavailable protein') |>
    mutate(correct = group == 'Deletion found' | group == 'Correct protein') |>
    group_by(tool, gene, correct) |>
    summarize(n = sum(n)) |>
    ungroup() |>
    group_by(tool, gene) |>
    reframe(correct = correct, perc = round(100 * n / sum(n), 1)) |>
    ungroup() |>
    as.data.frame()

