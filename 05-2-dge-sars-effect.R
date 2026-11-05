# Load libraries and custom functions ----

## Differential expression
library(DESeq2)

## Data manipulation
library(data.table)
library(dplyr)
library(tidyr)

## Visualization
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

# Objectives ----

# - Assess the effect of SARS-CoV-2 exposure (inactivated or infectious
#   virus) on gene expression in samples without plasma (design ~ sars).
# - Export DESeq2-normalized expression per sample and mean expression per
#   condition.
# - Generate PCA, DE tables and volcano plots.
# - Compare inactivated vs infectious SARS-CoV-2 effects (scatter plot,
#   log2FC correlation and Venn diagrams).

# User-defined parameters ----

## Global variables ----
onlyhuman <- FALSE

## Paths ----
output_path <- file.path(project_path, 'output', '05-dge')

## Load data ----
inputs <- load_dge_inputs(project_path, mapper)
dds <- inputs$dds
metadata <- inputs$metadata
features_db <- inputs$features_db
rm(inputs)

## Load palettes ----
source("r-utils/r-palettes.R")

# main ----

## Reformat data ----

# if required, select only human genes

genes <- select_genes(dds, features_db, onlyhuman)
mat <- genes$mat
genestag <- genes$genestag
rm(genes)

sars_path <- create_output_dirs(
  file.path(output_path, paste0('sars-effect-analysis', genestag)),
  subdirs = c('de-tables', 'volcano-plot'))

## Subset samples ----
# Goal: keep only samples without plasma treatment.

metadata_subset <- metadata[plasma == 'no_plasma', ]

# Use character IDs: indexing by a factor uses its integer codes, not labels.
target_samples <- as.character(metadata_subset$id)

#SOURCE: https://erilu.github.io/bulk-rnaseq-analysis/
dds_subset <- DESeqDataSetFromMatrix(
  countData = mat[, target_samples],
  colData = metadata_subset[, .(batch, sars)],
  design = ~ sars)

# Check that samples are in the same order within metadata and matrix.
stopifnot(
  identical(colnames(mat[, target_samples]), target_samples)
)

rm(dds, mat); gc()

# remove absent factor levels
vars <- c("sars", "batch")
for (v in vars) dds_subset[[v]] <- droplevels(dds_subset[[v]])
rm(vars, v)

## de analysis, design ~ sars ----
dds_subset <- DESeq(
  dds_subset,
  test = "Wald")
resultsNames(dds_subset)

### PCAs, vst-transformed ----
# Goal: run PCA on variance-stabilized data, selecting the top 500 most
# variable genes on the same scale.

pca_res <- run_pca_vst(dds_subset)

save_pca_plots(
  pca_res, metadata_subset,
  groups = 'sars',
  titles = 'PCA on VST-transformed data for sars effect analysis',
  out_dir = sars_path,
  file_prefix = 'pca-prcomp_vst')

### PCAs, log2-CPM values ----
# Goal: select the top 500 most variable genes and run the PCA on the same
# scale (TMM-normalized log2-CPM), so that gene selection is not driven by
# the mean-variance relationship of raw counts.

pca_res <- run_pca_lcpm(counts(dds_subset))

save_pca_plots(
  pca_res, metadata_subset,
  groups = 'sars',
  titles = 'PCA on log2-CPM data for sars effect analysis',
  out_dir = sars_path,
  file_prefix = 'pca-prcomp_log2CPM')

rm(pca_res)

### expression tables ----
# Goal: save DESeq2-normalized expression per sample, and the mean
# expression per condition (with plasma and sars annotation).

norm_mat <- counts(dds_subset, normalized = TRUE)

df_exp <- data.table(
  GENEID = row.names(norm_mat),
  norm_mat)

df_exp <- merge(
  df_exp,
  features_db[, .(GENEID, SYMBOL)],
  by = 'GENEID',
  all.x = TRUE)

fwrite(
  df_exp,
  file = file.path(sars_path, 'sars-effect_gene_expression_per_sample.csv'))

# Change to long format and add the condition of each sample.
df_exp_long <- df_exp |>
  pivot_longer(
    cols = -c(GENEID, SYMBOL),
    names_to = 'id',
    values_to = 'expression') |>
  left_join(
    metadata[, .(id = as.character(id), condition)],
    by = 'id')

# Calculate the mean expression per condition.
average_expression <- df_exp_long |>
  group_by(condition, GENEID, SYMBOL) |>
  summarise(avg_expression = mean(expression, na.rm = TRUE), .groups = 'drop') |>
  as.data.table()

average_expression <- merge(
  average_expression,
  unique(metadata[, .(condition, plasma, sars)]),
  by = 'condition',
  all.x = TRUE)

# Show results for SARS-CoV-2 N gene.
average_expression[SYMBOL == 'N']

fwrite(
  average_expression,
  file = file.path(sars_path, 'sars-effect_mean_gene_expression_per_condition.csv'))

rm(df_exp, df_exp_long, average_expression, norm_mat)

## sars effect ----

comparisons <- list(
  c("sars", "inactivated_sars", "no_sars"),
  c("sars", "infectious_sars", "no_sars"),
  c("sars", "infectious_sars", "inactivated_sars")
)

### de tables ----
invisible(lapply(comparisons, function(x) {
  write.csv(
    reformat_de_table(
      dds_subset,
      x,
      LFCCUT = NULL,
      PCUT = NULL,
      annotation = features_db),
    file = file.path(
      sars_path,
      'de-tables',
      paste0('de_sars-effect_', x[[2]], '_vs_', x[[3]], '.csv')),
    row.names = FALSE)
}))

### volcano plots ----
invisible(lapply(comparisons, function(x) {
  p <- do_volcanoplot(
    reformat_de_table(dds_subset, x, annotation = features_db),
    x_cutoff = 0.25,
    main = paste('sars effect,', x[[2]], 'vs', x[[3]]))

  save_svg_plot(
    p,
    file.path(
      sars_path, 'volcano-plot',
      paste0('de_volcano-plot_sars-effect_', x[[2]], '_vs_', x[[3]], '.svg')),
    width = 8, height = 8)
}))

## inactivated vs infectious sars comparison ----
# Goal: compare the effect of inactivated and infectious SARS-CoV-2 (both
# against no_sars), classifying genes by significance (padj < 0.05) and
# direction of change in each comparison.

diff_cmp <- merge(
  read_de_table(
    file.path(sars_path, 'de-tables',
              paste0('de_sars-effect_', comparisons[[1]][2], '_vs_', comparisons[[1]][3], '.csv')),
    cols = c('id', 'SYMBOL', 'padj', 'log2FoldChange'),
    suffix = 'inac_vs_no_sars'),
  read_de_table(
    file.path(sars_path, 'de-tables',
              paste0('de_sars-effect_', comparisons[[2]][2], '_vs_', comparisons[[2]][3], '.csv')),
    cols = c('id', 'SYMBOL', 'padj', 'log2FoldChange'),
    suffix = 'inf_vs_no_sars'),
  by = c('id', 'SYMBOL'),
  all = TRUE)

# NA padj (independent filtering / outliers) are treated as non-significant.
diff_cmp[, `:=`(
  sig.inac = !is.na(padj.inac_vs_no_sars) & padj.inac_vs_no_sars < 0.05,
  sig.inf  = !is.na(padj.inf_vs_no_sars)  & padj.inf_vs_no_sars  < 0.05)]

diff_cmp[, analysis := fcase(
  !sig.inac & !sig.inf, 'non-significant',
  sig.inac & sig.inf &
    sign(log2FoldChange.inac_vs_no_sars) == sign(log2FoldChange.inf_vs_no_sars), 'shared-gene',
  sig.inac & sig.inf, 'divergent genes',
  sig.inac, 'inactivated sars specific',
  default = 'infectious sars specific')]

print(table(diff_cmp$analysis))

### scatter plot ----
col_plot <- c(
  'non-significant' = 'grey',
  'shared-gene' = 'firebrick',
  'divergent genes' = 'black',
  'inactivated sars specific' = 'green4',
  'infectious sars specific' = 'dodgerblue')

diff_cmp[, gene_label := fifelse(
  abs(log2FoldChange.inac_vs_no_sars) >= 2 | abs(log2FoldChange.inf_vs_no_sars) >= 2,
  SYMBOL, '')]

p <- ggplot(diff_cmp, aes(
  x = log2FoldChange.inac_vs_no_sars,
  y = log2FoldChange.inf_vs_no_sars,
  label = gene_label,
  col = analysis)) +
  geom_point(size = 0.5) +
  ggrepel::geom_text_repel(size = 2, max.overlaps = 30) +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = 'dotted') +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dotted') +
  scale_color_manual(values = col_plot) +
  theme_minimal() +
  theme(legend.position = 'bottom')

save_svg_plot(
  p,
  file.path(sars_path, 'volcano-plot', 'scatterplot_sars_effect-inac.inf_vs_no-sars.svg'),
  width = 8, height = 6)

rm(p)

### correlation ----
# Goal: correlate log2FC of both comparisons using only human genes.

diff_cmp_human <- diff_cmp[id %in% features_db[species == 'homo sapiens', GENEID]]

corr_val <- cor.test(
  diff_cmp_human$log2FoldChange.inac_vs_no_sars,
  diff_cmp_human$log2FoldChange.inf_vs_no_sars)
print(corr_val)

rm(diff_cmp_human, corr_val)

### venn diagrams ----
# Goal: overlap of significant genes between both comparisons at two padj
# thresholds.

for (pcut in c(0.05, 0.001)) {
  vennData <- list(
    inac_no_sars = diff_cmp[padj.inac_vs_no_sars < pcut, id],
    inf_no_sars  = diff_cmp[padj.inf_vs_no_sars  < pcut, id])

  save_svg_plot(
    gplots::venn(vennData),
    file.path(
      sars_path,
      paste0('sars_effect_inactivated-infections_vs_no_sars_venndiag_PVAL', pcut, '.svg')),
    width = 8, height = 8)
}

rm(pcut, vennData, diff_cmp)

# session info ----

rm(list=ls())

gc()

sessioninfo::session_info()
