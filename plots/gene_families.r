if (!exists('.IMPORTED_COMMON')) { source('common.r') }

assign_gene_family <- function(df) {
    mutate(df, family = case_when(
        grepl('^MUC', locus) ~ 'MUC',
        grepl('^KCNK', locus) ~ 'KCNK',
        grepl('^CYP2', locus) ~ 'CYP2',
        locus == 'LPA' ~ 'LPA',
        grepl('^CFH', locus) ~ 'CFH',
        grepl('^FCGR', locus) ~ 'FCGR',
        T ~ NA)) |>
    filter(!is.na(family))
}

local({
    eval_dir <- file.path(proj_dir, 'summaries/eval')
    ver <- 'v0.14.2-4'
    summaries_fam <- data.frame()
    for (loo in c(F, T)) {
        summaries_fam <- read.csv(sprintf('%s/327%s_illumina_%s.csv.gz',
                eval_dir, ifelse(loo, 'LOO', ''), ver),
            sep = '\t', comment = '#') |>
            mutate(tech = ifelse(loo, 'Locityper LOO', 'Locityper')) |>
            assign_gene_family() %>%
            smartbind(summaries_fam, .)
    }
    
    summaries_fam <- read.csv(
        sprintf('%s/comparison/NYGC/summaries/nygc.eval.csv.gz', proj_dir),
            sep = '\t', comment = '#') |>
        mutate(tech = '1KGP') |>
        assign_gene_family() %>%
        smartbind(summaries_fam, .)
    summaries_fam <<- preproc_summary(summaries_fam)
})

MAX_QV <- 50
colors <- ggthemes::tableau_color_pal('Classic 10 Medium')(10)[c(5, 2, 3, 4, 1)]

aggr <- mutate(summaries_fam, qv = pmin(qv, MAX_QV)) %>%
    aggregate(qv ~ locus + tech + family, ., mean) |>
    pivot_wider(names_from = tech, values_from = qv)

ggplot(aggr) +
    geom_abline() +
    geom_point(aes(`1KGP`, `Locityper LOO`, color = family)) +
    scale_color_manual(values = c(colors, 'black')) +
    theme_bw()


#########

library(bignum)
geom_mean <- function(x) { exp(mean(log(as_bigfloat(x)))) }
ext_geom_mean <- function(x, delta) { exp(mean(log(as_bigfloat(x) + delta))) - delta }
est_delta <- function(x, eps) {
    x <- bigfloat(x)
    pos_x <- x[x > 0]
    if (length(unique(pos_x)) <= 1) {
        return(NA)
    }

    TEN <- bigfloat(10)
    inner <- function(log_delta) {
        ((1 + eps) * geom_mean(pos_x) - ext_geom_mean(pos_x, TEN^log_delta)) |>
            as.numeric() |>
            suppressWarnings()
    }
    LEFT <- -20
    RIGHT <- 10
    if (inner(RIGHT) > 0) {
        return(10 ^ RIGHT)
    } else if (inner(LEFT) < 0) {
        return(10 ^ LEFT)
    }

    10 ^ uniroot(inner, c(LEFT, RIGHT), tol = 1e-10)$root
}

colors <- ggthemes::tableau_color_pal('Classic 10 Medium')(10)[c(5, 2, 3, 4, 1)]

EPS <- 1e-5
deltas <- aggregate(div ~ locus + tech + family, summaries_fam, est_delta, EPS)
locus_deltas <- aggregate(div ~ locus, summaries_fam, est_delta, 1e-4) |>
    rename(delta = div)

aggr <- aggregate(qv ~ locus + tech + family, summaries_fam, function(x) mean(pmin(x, 160))) |>
    mutate(div = NA)

locus_deltas
(delta4 <- est_delta(summaries_fam$div, 1e-4))
(delta5 <- est_delta(summaries_fam$div, 1e-5))
aggr <- aggregate(div ~ locus + tech + family, summaries_fam, ext_geom_mean, delta4) |>
    mutate(
        div = pmax(0, suppressWarnings(as.numeric(div))),
        qv = pmin(160, -10 * log10(div)))

aggr <- left_join(summaries_fam, locus_deltas, by = 'locus') |>
    group_by(locus, tech, family) |>
    summarize(
        div = ext_geom_mean(div, unique(delta)),
        .groups = 'keep'
    ) |>
    ungroup() |>
    mutate(
        div = pmax(0, suppressWarnings(as.numeric(div))),
        qv = pmin(60, -10 * log10(div))
    )

pivot_wider(aggr, id_cols = -div, names_from = tech, values_from = qv) |>
ggplot() +
    geom_abline() +
    geom_point(aes(`1KGP`, `Locityper LOO`, color = family)) +
    scale_color_manual(values = c(colors, 'black')) +
    theme_bw()
filter(aggr, locus == 'CYP2J2')
filter(summaries_fam, locus == 'CYP2J2') |> select(tech, div, qv) |> filter(div != 0)


