# Load libraries ----

## Differential expression
library(DESeq2)

## Data manipulation
library(data.table)

## Visualization
library(ggplot2)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")

# Objectives ----

# - Plot the normalized counts of the top variable genes (vst data) by
#   SARS-CoV-2 and plasma condition.
# - Check a housekeeping gene (ACTB) and the SARS-CoV-2 transcript counts per
#   condition.

# User-defined parameters ----

## Global variables ----

N_TOP_GENES  <- 5      # top variable genes to plot
HOUSEKEEPING <- 'ACTB' # housekeeping gene to check

## Paths ----

output_path <- file.path(
  project_path,
  'output',
  '07-overview',
  'sars-infection')
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## Load data ----

inputs <- load_dge_inputs(project_path, mapper)
dds <- inputs$dds
metadata <- inputs$metadata
features_db <- inputs$features_db
setkey(features_db, 'GENEID')
rm(inputs)

## Load palettes ----

source("r-utils/r-palettes.R")

## Functions ----

# Sample metadata with id as character (for joins with matrix column names).
sample_info <- function(cols) {
  md <- metadata[, c('id', cols), with = FALSE]
  md[, id := as.character(id)]
  md
}

# Counts (ComBat-seq adjusted) of a set of genes, one row per sample,
# matched to the metadata by sample id.
counts_per_sample <- function(dds, gene_ids) {
  mat <- floor(assay(dds)[gene_ids, , drop = FALSE])

  dt <- as.data.table(t(mat))
  setnames(dt, make.unique(features_db[.(gene_ids), GENELAB]))
  dt[, id := colnames(mat)]

  dt <- merge(dt, sample_info('conditions'), by = 'id', sort = FALSE)
  setcolorder(dt, 'id')
  dt
}

plot_counts_by_condition <- function(dt, y) {
  ggplot(dt, aes(x = conditions, y = .data[[y]], color = conditions)) +
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Main processing ----

## Variance stabilization ----
# Goal: vst-transform the data and select the most variable genes.
# The counts are already batch-corrected with ComBat-seq (03-batch-correction),
# so no further batch correction is applied. Removing the batch effect again
# with limma::removeBatchEffect() gave the same top variable genes.

design(dds) <- ~ plasma + sars

vsd <- DESeq2::vst(
  dds,           # normalized counts cannot be used as input
  blind = FALSE, # use the design to estimate the dispersion trend
  fitType = 'parametric',
  nsub = 1000)

## Top variable genes ----

rv <- rowVars(assay(vsd))
target_ids <- rownames(vsd)[order(rv, decreasing = TRUE)[seq_len(N_TOP_GENES)]]

top_genes <- data.table(
  rank = seq_len(N_TOP_GENES),
  GENEID = target_ids,
  GENELAB = features_db[.(target_ids), GENELAB])
print(top_genes)

fwrite(
  top_genes,
  file.path(output_path, paste0('top', N_TOP_GENES, '_variable_genes.tsv')),
  sep = '\t')

rm(vsd, rv)

## Count plots of top variable genes ----

count_plot_path <- file.path(output_path, 'count-plot')
dir.create(count_plot_path, recursive = TRUE, showWarnings = FALSE)

for (gene in target_ids) {
  
  d <- plotCounts(
    dds,
    gene = gene,
    intgroup = 'plasma',
    returnData = TRUE) |>
    as.data.table(keep.rownames = 'id')
  d <- merge(d, sample_info(c('conditions', 'sars')), by = 'id')
  
  gene_name <- features_db[.(gene), GENELAB]
  
  p <- ggplot(d, aes(x = sars, y = count, color = plasma, group = plasma)) +
    labs(
      x = 'Condition',
      y = 'Normalized counts',
      title = paste0(gene_name, ' is in top', N_TOP_GENES, ' variable gene')) +
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
    stat_summary(aes(y = count), fun = mean, geom = 'line') +
    scale_y_log10() +
    scale_color_manual(values = plasma_palette) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 22),
      strip.text = element_text(size = 15))
  # SOURCE: https://bookdown.org/cmecklin/sta265notes/anova-with-interaction.html
  
  save_svg_plot(
    p,
    file.path(count_plot_path, paste0('countplot_', gene_name, '.svg')),
    width = 7, height = 7)
}

## Housekeeping gene counts ----
# Goal: check that the housekeeping gene is stable across conditions.

hk_id <- features_db[SYMBOL == HOUSEKEEPING, GENEID]
dt_hk <- counts_per_sample(dds, hk_id)

fwrite(
  dt_hk,
  file = file.path(output_path, paste0(tolower(HOUSEKEEPING), '_genes_count.filter.combat.tsv')),
  sep = '\t')

save_svg_plot(
  plot_counts_by_condition(dt_hk, features_db[.(hk_id[1]), GENELAB]),
  file.path(output_path, paste0(tolower(HOUSEKEEPING), '_count.filter.combat.svg')),
  width = 8, height = 6)

## SARS-CoV-2 transcript counts ----
# Goal: count the viral reads per gene and in total for every sample.

sars_ids   <- features_db[species == 'sars-cov-2', GENEID]
sars_genes <- make.unique(features_db[.(sars_ids), GENELAB])

dt_infection <- counts_per_sample(dds, sars_ids)
dt_infection[, total := rowSums(.SD), .SDcols = sars_genes]
setcolorder(dt_infection, c('id', sars_genes, 'total', 'conditions'))

fwrite(
  dt_infection,
  file = file.path(output_path, 'sars_genes_count.filter.combat.tsv'),
  sep = '\t')

save_svg_plot(
  plot_counts_by_condition(dt_infection, 'total'),
  file.path(output_path, 'sars_count.filter.combat.svg'),
  width = 8, height = 6)

dt_infection_long <- melt(
  dt_infection,
  id.vars = c('id', 'conditions'),
  variable.name = 'genes',
  value.name = 'Counts',
  na.rm = TRUE)

p <- plot_counts_by_condition(dt_infection_long, 'Counts') +
  facet_wrap(~ genes, ncol = 2)

save_svg_plot(
  p,
  file.path(output_path, 'sars_count-per-gene.filter.combat.svg'),
  width = 7, height = 10)

# Session info ----

rm(list = ls())
gc()
sessioninfo::session_info()
