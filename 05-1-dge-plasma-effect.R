# Load libraries and custom functions ----

## Differential expression
library(DESeq2)

## Data manipulation
library(data.table)
library(dplyr)
library(tidyr)

## Visualization
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(RColorBrewer)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

# Objectives ----

# - Assess the effect of plasma (COVID-19 vs healthy donors) on gene
#   expression, adjusting for SARS-CoV-2 exposure (design ~ sars + plasma).
# - Export DESeq2-normalized expression per sample and mean expression per
#   condition.
# - Generate PCA, DE tables and volcano plots.

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

plasma_path <- create_output_dirs(
  file.path(output_path, paste0('plasma-effect-analysis', genestag)),
  subdirs = c('de-tables', 'volcano-plot'))

## Subset samples ----
# Goal: keep only samples treated with plasma (healthy or COVID-19 donors).

metadata_subset <- metadata[plasma != 'no_plasma', ]

# Use character IDs: indexing by a factor uses its integer codes, not labels.
target_samples <- as.character(metadata_subset$id)

#SOURCE: https://erilu.github.io/bulk-rnaseq-analysis/
dds_subset <- DESeqDataSetFromMatrix(
  countData = mat[, target_samples],
  colData = metadata_subset[, .(batch, sars, plasma)],
  design = ~ sars + plasma)

# Check that samples are in the same order within metadata and matrix.
stopifnot(
  identical(colnames(mat[, target_samples]), target_samples)
)

rm(dds, mat); gc()

# remove absent factor levels
vars <- c("sars", "plasma", "batch")
for (v in vars) dds_subset[[v]] <- droplevels(dds_subset[[v]])
remove(vars, v)

## de analysis, design ~ sars + plasma ----
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
  groups = c('plasma', 'sars'),
  titles = c(
    plasma = 'PCA on VST-transformed data for plasma effect analysis',
    sars   = 'PCA on VST-transformed data for sars effect analysis'),
  out_dir = plasma_path,
  file_prefix = 'pca-prcomp_vst')

### PCAs, log2-CPM values ----
# Goal: select the top 500 most variable genes and run the PCA on the same
# scale (TMM-normalized log2-CPM), so that gene selection is not driven by
# the mean-variance relationship of raw counts.

pca_res <- run_pca_lcpm(counts(dds_subset))

save_pca_plots(
  pca_res, metadata_subset,
  groups = c('plasma', 'sars'),
  titles = c(
    plasma = 'PCA on log2-CPM data for plasma effect analysis',
    sars   = 'PCA on log2-CPM data for sars effect analysis'),
  out_dir = plasma_path,
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
  file = file.path(plasma_path, 'plasma-effect_gene_expression_per_sample.csv'))

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
  file = file.path(plasma_path, 'plasma-effect_mean_gene_expression_per_condition.csv'))

rm(df_exp, df_exp_long, average_expression, norm_mat)

## plasma effect ----

comparisons <- list(
  c("plasma", "plasma_covid", "plasma_healthy")
)

res <- results(dds_subset, contrast = comparisons[[1]])
table(res$padj < 0.05)

# NAs: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
res[is.na(res$padj), ]
rm(res)

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
      plasma_path,
      'de-tables',
      paste0('de_plasma-effect_', x[[2]], '_vs_', x[[3]], '.csv')),
    row.names = FALSE)
}))

### volcano plots ----
invisible(lapply(comparisons, function(x) {
  p <- do_volcanoplot(
    reformat_de_table(dds_subset, x, annotation = features_db),
    x_cutoff = 0.25,
    main = paste('plasma effect,', x[[2]], 'vs', x[[3]]))

  save_svg_plot(
    p,
    file.path(
      plasma_path, 'volcano-plot',
      paste0('de_volcano-plot_plasma-effect_', x[[2]], '_vs_', x[[3]], '.svg')),
    width = 8, height = 8)
}))

# session info ----

rm(list=ls())

gc()

sessioninfo::session_info()
