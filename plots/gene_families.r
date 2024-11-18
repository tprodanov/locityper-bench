if (!exists('.IMPORTED_COMMON')) { source('common.r') }

assign_gene_family <- function(df) {
    mutate(df, family = case_when(
        grepl('^MUC', locus) ~ 'MUC',
        grepl('^CYP2', locus) ~ 'CYP2',
        grepl('^CFH', locus) ~ 'CFH',
        grepl('^FCGR', locus) ~ 'FCGR',
        T ~ NA)) |>
    filter(!is.na(family))
}

local({
    ver <- 'v0.17.0'
    summaries_fam <- data.frame()
    for (loo in c(F, T)) {
        summaries_fam <- read.csv(sprintf('%s/%s_illumina_%s.csv.gz',
                eval_dir, ifelse(loo, 'loo', 'full'), ver),
            sep = '\t', comment = '#') |>
            mutate(tech = ifelse(loo, 'Locityper LOO', 'Locityper')) |>
            assign_gene_family() %>%
            smartbind(summaries_fam, .)
    }

    summaries_fam <- read.csv(
        file.path(proj_dir, 'comparison/NYGC/summaries/nygc.eval.csv.gz'),
            sep = '\t', comment = '#') |>
        mutate(tech = '1KGP') |>
        assign_gene_family() %>%
        smartbind(summaries_fam, .)
    thresholds <- c(0, 17, 23, 33)
    summaries_fam <<- group_qvs(summaries_fam, thresholds) |>
        bound_qv() |>
        assign_filters()
})
samples_39 <- readLines(file.path(proj_dir, 'comparison/NYGC/samples_39.txt'))
summaries_fam <- filter(summaries_fam, sample %in% samples_39)

#########

scales::show_col(classic10 <- ggthemes::tableau_color_pal('Classic 10 Medium')(10))
colors <- c(
    'CFH' = classic10[2],
    'CYP2' = classic10[9],
    'FCGR' = classic10[4],
    'MUC' = classic10[1]
)

aggr <- aggregate(qv ~ tech + locus + family, summaries_fam, mean)
waggr <- pivot_wider(aggr, names_from = 'tech', values_from = 'qv') |>
    rename(nygc = `1KGP`, full = Locityper, loo = `Locityper LOO`) |>
    mutate(improv_full = full - nygc, improv_loo = loo - nygc) |>
    arrange(improv_loo)

annot_y <- max(waggr$loo)
lines <- c(10, 20, 30)
ggplot(waggr) +
    geom_abline(intercept = lines, color = 'gray90') +
    annotate('label', x = annot_y - lines, y = annot_y, label = lines,
        family = 'Overpass',
        color = 'gray30', size = 3,
        label.size = 0, label.padding = unit(0.15, "lines")) +
    geom_abline(color = 'gray50') +

    geom_point(aes(nygc, loo, fill = family),
        alpha = 0.9, shape = 21, color = 'gray30', stroke = 0.5, size = 2.8) +
    coord_fixed() +
    scale_x_continuous('1KGP accuracy (QV)',
        breaks = seq(0, 100, 10), expand = expansion(add = 1.05)) +
    scale_y_continuous('Locityper accuracy (QV)',
        breaks = seq(0, 100, 10), expand = expansion(add = 1.05)) +
    scale_fill_manual(NULL, values = colors) +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.9, 0.0),
        legend.text = element_text(margin = margin(l = -2, r = 0)),
        legend.margin = margin(b = 4, l = 1, r = 1, t = -2),
        legend.justification = 'bottom',
        plot.margin = margin(l = 1),
        axis.title.y = element_text(hjust = 0),
    )
ggsave(file.path(plots_dir, 'gene_families.png'),
    width = 5.3, height = 4, dpi = 600, scale = 0.63)

#############

secreted <- sprintf('MUC%s', c('6', '2', '19', '5AC', '5B', '7'))
mucs <- filter(waggr, family == 'MUC') |>
    arrange(desc(improv_loo)) |>
    mutate(
        type = factor(ifelse(locus %in% secreted, 'secreted', 'tethered'), levels = c('tethered', 'secreted')),
        locus = factor(locus, levels = unique(locus))) |>
    select(locus, type, improv_full, improv_loo) |>
    pivot_longer(c(improv_full, improv_loo), names_to = 'database', values_to = 'improv')  |>
    mutate(database = recode_factor(database,
        'improv_loo' = 'leave-one-out', 'improv_full' = 'full database'))

ggplot(mucs) +
    geom_bar(aes(locus, improv, group = database), fill = 'white',
        stat = 'identity', width = 0.8,
        position = position_dodge2(padding = 0)) +
    geom_bar(aes(locus, improv, fill = type, group = database, alpha = database),
        stat = 'identity', width = 0.8,
        position = position_dodge2(padding = 0)) +
    scale_x_discrete(NULL) +
    scale_y_continuous('Accuracy improvement<br/>(QV<sub>Locityper</sub> – QV<sub>1KGP</sub>)',
        expand = expansion(add = 1)) +
    scale_fill_manual(NULL, values = c(rgb(109, 97, 168), '#1ca449')) +
    scale_alpha_manual(NULL, values = c(1, 0.5)) +
    guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
    theme_bw() +
    theme(
        text = element_text(family = 'Overpass'),
        panel.border = element_blank(),
        axis.text.x = element_text(
            size = 9, angle = 35, vjust = 1, hjust = 1, face = 'italic'),
        axis.title.y = ggtext::element_markdown(
            size = 11, lineheight = 1.2, margin = margin(t = 10, r = 2)),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.spacing.y = unit(-0.1, 'lines'),
        legend.text = element_text(margin = margin(t = 1, b = 1, r = 20, l = 5)),
        legend.key.width = unit(1, 'lines'),
        legend.key.height = unit(0.4, 'lines'),
        legend.position = 'inside',
        legend.position.inside = c(0.91, 0.735),
        plot.margin = margin(l = 1),
    )
ggsave(file.path(plots_dir, 'improv_MUC.png'),
    width = 11, height = 5, dpi = 600, scale = 0.48)

#####################

aggregate(improv_loo ~ family, waggr, mean)

filter(waggr, family == 'MUC')
filter(waggr, family == 'FCGR')
filter(waggr, family == 'CFH')
filter(waggr, family == 'CYP2')
