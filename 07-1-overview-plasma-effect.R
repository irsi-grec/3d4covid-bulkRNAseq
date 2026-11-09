# Load libraries ----

## Annotation
library(org.Hs.eg.db)

## Data manipulation
library(data.table)

## Visualization
library(ggplot2)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-plots.R")

# Objectives ----

# - Plot the log2FC of the DE genes annotated to selected GO terms
#   (plasma_covid vs plasma_healthy).
# - Plot the mean expression per condition of the genes in those GO terms.
#
# Where two options exist, both are run and labelled:
# - [SUBMITTED] reproduces the version sent for review.
# - [REVISED]   corrected/alternative version.
#
# Options (all combinations are run):
# - GO keytype:
#   - [SUBMITTED] 'GO': genes annotated directly to the GO term.
#   - [REVISED]   'GOALL': also genes annotated to its descendant terms.
# - Gene matching between GO genes (org.Hs.eg.db) and DE tables (EnsDb v86):
#   - [SUBMITTED] 'SYMBOL': symbols may differ between both annotations.
#   - [REVISED]   'ENSEMBL': Ensembl gene IDs.
# The [SUBMITTED] combination (GO + SYMBOL) keeps the original file names;
# the other ones add the suffix '_<keytype>-<matching>'.

# User-defined parameters ----

## Global variables ----

# Comparisons to plot (one heatmap column per comparison).
TARGET_CMP <- c('plasma_covid_vs_plasma_healthy')
# c('inactivated_sars_vs_no_sars', 'infectious_sars_vs_no_sars', 'infectious_sars_vs_inactivated_sars')

# GO category to plot (column 'categorie' in GO_list.csv).
GO_CATEGORY <- 'Vasculature'

PADJ_CUT  <- 0.05
N_COLUMNS <- 4 # number of columns in which the heatmap is split

# Options to run (see Objectives). The first combination is [SUBMITTED].
GO_OPTIONS <- CJ(
  keytype  = c('GO', 'GOALL'),
  match_by = c('SYMBOL', 'ENSEMBL'),
  sorted = FALSE)
GO_OPTIONS[, suffix := fifelse(
  keytype == 'GO' & match_by == 'SYMBOL', '', paste0('_', keytype, '-', match_by))]

## Paths ----

data_path <- file.path(
  project_path,
  'output',
  '01-data-formating')

de_path <- file.path(
  project_path,
  'output',
  '05-dge',
  'plasma-effect-analysis')

go_file <- file.path(
  project_path,
  'data',
  'GO',
  'GO_list.csv')

output_path <- file.path(
  project_path,
  'output',
  '07-overview',
  'plasma-effect',
  'selected-GO')
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## Genes database ----

features_db <- readRDS(file.path(data_path, 'genesDb.rds'))
setkey(features_db, 'GENEID')

## Functions ----

# Assign n genes to m heatmap columns, filling the first columns first
# (e.g. 10 genes in 4 columns -> 3, 3, 2, 2).
div_and_list <- function(n, m) {
  sort(rep_len(seq_len(m), n))
}

# Add the heatmap column of every gene (alphabetical order) and set SYMBOL
# as an ordered factor.
add_heatmap_columns <- function(dt, m) {
  symbols  <- sort(unique(as.character(dt$SYMBOL)))
  gene2col <- setNames(div_and_list(length(symbols), m), symbols)

  dt[, ncolumn := gene2col[as.character(SYMBOL)]]
  dt[, SYMBOL := factor(SYMBOL, levels = symbols)]
  dt
}

# Tile plot for the genes of one heatmap column.
plot_tile_column <- function(dt, xcolumn, x, fill, fill_range, label = NULL) {
  p <- ggplot(dt[ncolumn == xcolumn], aes(x = .data[[x]], y = SYMBOL, fill = .data[[fill]])) +
    geom_tile() +
    scale_y_discrete(limits = rev) + # A at the top
    scale_fill_distiller(palette = 'RdBu', limits = c(-fill_range, fill_range)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
  
  if (!is.null(label)) {
    p <- p + geom_text(aes(label = .data[[label]]), color = 'black', size = 4, vjust = 0.75)
  }
  p
}

# Arrange the heatmap columns side by side with a common legend.
plot_tile_grid <- function(dt, m, ...) {
  plots_list <- lapply(seq_len(m), function(i) plot_tile_column(dt, xcolumn = i, ...))
  ggpubr::ggarrange(
    plotlist = plots_list,
    ncol = m,
    common.legend = TRUE,
    legend = 'right')
}

# Genes annotated to the target GO terms (columns: GO, ENSEMBL, SYMBOL).
get_go_genes <- function(go_ids, keytype) {
  dt <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_ids,
    columns = c('ENSEMBL', 'SYMBOL'),
    keytype = keytype) |>
    as.data.table()

  setnames(dt, keytype, 'GO')
  unique(dt[GO %in% go_ids, .(GO, ENSEMBL, SYMBOL)])
}

read_de_table_cmp <- function(cmp_id) {
  de_file <- file.path(de_path, 'de-tables', paste0('de_plasma-effect_', cmp_id, '.csv'))
  if (!file.exists(de_file)) {
    warning('DE table not found: ', de_file)
    return(NULL)
  }
  dt <- fread(de_file)
  dt <- dt[!is.na(padj) & !is.na(SYMBOL)]
  dt[, cmp := cmp_id]
  dt
}

# Main processing ----

## Load data ----

go_list   <- fread(go_file)
go_target <- go_list[categorie == GO_CATEGORY, unique(goid)]

de_table_all <- rbindlist(lapply(TARGET_CMP, read_de_table_cmp))
de_table_all[, cmp := factor(cmp, levels = TARGET_CMP)]
de_table_all[, sig := fifelse(padj < PADJ_CUT, '*', '')]

mean_file <- file.path(de_path, 'plasma-effect_mean_gene_expression_per_condition.csv')
mean_exp_all <- if (file.exists(mean_file)) fread(mean_file)[!is.na(SYMBOL)] else NULL

## Plots per GO option ----

option_summary <- list()

for (i in seq_len(nrow(GO_OPTIONS))) {

  opt <- GO_OPTIONS[i]
  message('GO keytype: ', opt$keytype, ' | matching by: ', opt$match_by)

  ### Target GO genes ----

  go2genes <- get_go_genes(go_target, opt$keytype)

  if (opt$match_by == 'SYMBOL') {
    # [SUBMITTED] symbols from org.Hs.eg.db.
    in_go_de   <- function(dt) dt$SYMBOL %in% go2genes[!is.na(SYMBOL), SYMBOL]
    go_symbols <- unique(go2genes[!is.na(SYMBOL), SYMBOL])
  } else {
    # [REVISED] Ensembl IDs; symbols taken from EnsDb v86 (as in the DE tables).
    go_ids     <- unique(go2genes[!is.na(ENSEMBL), ENSEMBL])
    in_go_de   <- function(dt) dt$id %in% go_ids
    go_symbols <- features_db[.(go_ids), SYMBOL, nomatch = NULL] |> na.omit() |> unique()
  }

  ### DE genes in target GOs ----
  # Goal: plot the log2FC of the genes in the target GO terms that are
  # significant in at least one comparison.

  de_table <- de_table_all[in_go_de(de_table_all)]
  de_table <- de_table[SYMBOL %in% de_table[padj < PADJ_CUT, SYMBOL]]
  de_table <- add_heatmap_columns(de_table, N_COLUMNS)

  option_summary[[i]] <- data.table(
    keytype = opt$keytype,
    match_by = opt$match_by,
    go_genes = uniqueN(go2genes$ENSEMBL),
    de_genes_plotted = uniqueN(de_table$SYMBOL))

  if (nrow(de_table) > 0) {
    p <- plot_tile_grid(
      de_table,
      m = N_COLUMNS,
      x = 'cmp',
      fill = 'log2FoldChange',
      fill_range = max(abs(de_table$log2FoldChange)),
      label = 'sig')

    file_base <- file.path(output_path, paste0('go_plasma-effect_log2FC_GO-vascular', opt$suffix))
    save_svg_plot(p, paste0(file_base, '.svg'), width = 8, height = 10)
    fwrite(de_table, paste0(file_base, '.tsv'), sep = '\t')
  }

  ### Mean expression in target GOs ----
  # Goal: plot the mean expression per condition of all genes in the target
  # GO terms.

  if (!is.null(mean_exp_all)) {
    mean_exp <- mean_exp_all[SYMBOL %in% go_symbols]
    mean_exp[, sars := factor(sars, levels = c('no_sars', 'inactivated_sars', 'infectious_sars'))]
    mean_exp <- add_heatmap_columns(mean_exp, N_COLUMNS)

    p <- plot_tile_grid(
      mean_exp,
      m = N_COLUMNS,
      x = 'condition',
      fill = 'avg_expression',
      fill_range = max(abs(mean_exp$avg_expression)))

    save_svg_plot(
      p,
      file.path(output_path, paste0('go_plasma-effect_mean-gene-exp_GO-vascular', opt$suffix, '.svg')),
      width = 12, height = 40)
  }
}

## Summary of GO options ----

option_summary <- rbindlist(option_summary)
print(option_summary)

fwrite(
  option_summary,
  file.path(output_path, 'go_plasma-effect_GO-vascular_options-summary.tsv'),
  sep = '\t')

# Session info ----

rm(list = ls())
gc()
sessioninfo::session_info()
