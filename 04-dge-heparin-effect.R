# Load libraries and custom functions ----

## Differential expression
library(DESeq2)
library(gprofiler2)

library(data.table)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

library(agricolae)

## Data manipulation
library(data.table)
library(tidyr)

## Visualization
library(ggplot2)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

set.seed(123)

# User-defined parameters ----

## Global variables ----
onlyhuman <- FALSE

## Paths ----
output_path <- file.path(project_path, 'output', '04-dge-heparin-effect')

## Load data ----
inputs <- load_dge_inputs(project_path, mapper)
dds <- inputs$dds
metadata <- inputs$metadata
features_db <- inputs$features_db
rm(inputs)

## Load palettes ----
source("r-utils/r-palettes.R")

## Functions ----

compare_models <- function(x, padj_threshold = NULL, log2fc_threshold = NULL) {
  
  col_plot <- c(
    'shared-gene' = 'firebrick',
    'divergent genes' = 'grey',
    'simple model specific' = 'green4',
    'additive model specific' = 'dodgerblue')
  
  cmp_id   <- paste0(x[[2]], '_vs_', x[[3]])
  out_base <- file.path(output_path, 'model-comparison',
                        paste0('sars_effect-simple_vs_additive_model-', cmp_id))
  sfx <- make_suffix(padj_threshold, log2fc_threshold)
  
  sim.mod <- read_de_table(
    file.path(output_path, 'de-tables',
              paste0('de_heparin-effect_simple-model_', cmp_id, '.csv')),
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    cols = c('id', 'SYMBOL', 'log2FoldChange'),
    suffix = 'simplemod')
  
  add.mod <- read_de_table(
    file.path(output_path, 'de-tables',
              paste0('de_heparin-effect_additive-model_', cmp_id, '.csv')),
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    cols = c('id', 'SYMBOL', 'log2FoldChange'),
    suffix = 'additivemod')
  
  # Do Venn diagram
  vennData <- list(
    simple.model   = sim.mod$id,
    additive.model = add.mod$id)
  
  save_svg_plot(
    gplots::venn(vennData),
    paste0(out_base, '_venndiag', sfx, '.svg'), width = 8, height = 8)
  
  # LFC comparison
  diff.cmp <- merge(sim.mod, add.mod, by = c('id', 'SYMBOL'), all = TRUE)
  setnafill(diff.cmp, fill = 0, cols = c('log2FoldChange.simplemod', 'log2FoldChange.additivemod'))
  
  diff.cmp[, analysis := fcase(
    log2FoldChange.simplemod == 0 & log2FoldChange.additivemod == 0, 'both zero',
    log2FoldChange.simplemod * log2FoldChange.additivemod > 0,       'shared-gene',
    log2FoldChange.simplemod != 0 & log2FoldChange.additivemod == 0, 'simple model specific',
    log2FoldChange.simplemod == 0 & log2FoldChange.additivemod != 0, 'additive model specific',
    default = 'divergent genes')]
  print(table(diff.cmp$analysis))
  
  fwrite(diff.cmp, paste0(out_base, sfx, '.csv'))
  
  diff.cmp <- diff.cmp[analysis != 'both zero']
  print(cor.test(diff.cmp$log2FoldChange.simplemod, diff.cmp$log2FoldChange.additivemod))
  
  diff.cmp[, gene_label := fifelse(
    abs(log2FoldChange.simplemod) >= 2 | abs(log2FoldChange.additivemod) >= 2, SYMBOL, '')]
  
  p <- ggplot(diff.cmp, aes(
    x = log2FoldChange.simplemod, y = log2FoldChange.additivemod,
    label = gene_label, col = analysis)) +
    geom_point(size = 0.5) +
    ggrepel::geom_text_repel(size = 3) +
    scale_color_manual(values = col_plot) +
    labs(title = paste(x[[2]], 'vs', x[[3]])) +
    theme_minimal() +
    theme(legend.position = 'bottom')
  
  save_svg_plot(p, paste0(out_base, sfx, '.svg'), width = 8, height = 6) 
  
  invisible(diff.cmp)
}

# main ----

## reformat data ----

# Create a DGEList object.
# if required, select only human genes

genes <- select_genes(dds, features_db, onlyhuman)
mat <- genes$mat
genestag <- genes$genestag
rm(genes)
if (onlyhuman) output_path <- paste0(output_path, genestag)

create_output_dirs(
  output_path,
  subdirs = c('model-comparison', 'de-tables', 'volcano-plot'))

# remove absent factor levels
vars <- c("heparin", "sars", "plasma", "batch")
for (v in vars) dds[[v]] <- droplevels(dds[[v]])
rm(vars)

## subset samples ----
target_samples <- metadata[plasma == 'no_plasma', id]
metadata_subset <- metadata[plasma == 'no_plasma',]

#SOURCE: https://erilu.github.io/bulk-rnaseq-analysis/
dds_subset <- DESeqDataSetFromMatrix(
  countData = mat[, target_samples],
  colData = metadata_subset[, .(heparin, sars)],
  design = ~ heparin + sars)

## de analysis, design ~  sars ----
design(dds_subset) <- ~  sars
dds_subset <- DESeq(dds_subset)
resultsNames(dds_subset)

comparisons <- list(
  c("sars", "inactivated_sars", "no_sars"),  
  c("sars", "infectious_sars", "no_sars"),
  c("sars", "infectious_sars", "inactivated_sars")
)

invisible(lapply(comparisons, function(x) {
  write.csv(
  reformat_de_table(
    dds_subset,
    x,
    LFCCUT = 0.25,
    PCUT = 0.05,
    annotation = features_db),
  file = file.path(
    output_path,
    'de-tables', 
    paste0('de_heparin-effect_simple-model_', x[[2]], '_vs_', x[[3]], '_LFCUT0.25_PCUT0.05.csv')),
  row.names = FALSE)
}))

invisible(lapply(comparisons, function(x) {
  write.csv(
  reformat_de_table(
    dds_subset,
    x,
    LFCCUT = NULL,
    PCUT = NULL,
    annotation = features_db),
  file = file.path(
    output_path,
    'de-tables',
    paste0('de_heparin-effect_simple-model_', x[[2]], '_vs_', x[[3]], '.csv')),
  row.names = FALSE)
}))

### volcano plots ----
invisible(lapply(comparisons, function(x) {
  p <- do_volcanoplot(
    reformat_de_table(dds_subset, x, annotation = features_db),
    x_cutoff = 0.25,
    main = paste(x[[2]], 'vs', x[[3]]))
  
  save_svg_plot(
    p,
    file.path(
      output_path, 'volcano-plot',
      paste0('de-volcano-plot_heparin-effect_simple-model_', x[[2]], '_vs_', x[[3]], '.svg')),
    width = 8, height = 8)
}))

## de analysis, design ~  sars + heparin ----
design(dds_subset) <- ~ sars + heparin
dds_subset <- DESeq(dds_subset)
resultsNames(dds_subset)

comparisons <- list(
  c("heparin", "heparin", "no_heparin"),
  c("sars", "inactivated_sars", "no_sars"),
  c("sars", "infectious_sars", "no_sars"),
  c("sars", "infectious_sars", "inactivated_sars")
)

invisible(lapply(comparisons, function(x) {
  write.csv(
  reformat_de_table(
    dds_subset,
    x,
    LFCCUT = 0.25,
    PCUT = 0.05,
    annotation = features_db),
  file = file.path(output_path, 'de-tables', paste0('de_heparin-effect_additive-model_', x[[2]], '_vs_', x[[3]], '_LFCUT0.25_PCUT0.05.csv')),
  row.names = FALSE)
}))

invisible(lapply(comparisons, function(x) {
  write.csv(
  reformat_de_table(
    dds_subset,
    x,
    LFCCUT = NULL,
    PCUT = NULL,
    annotation = features_db),
  file = file.path(output_path, 'de-tables', paste0('de_heparin-effect_additive-model_', x[[2]], '_vs_', x[[3]], '.csv')),
  row.names = FALSE)
}))

### volcano plots ----

invisible(lapply(comparisons, function(x) {
  p <- do_volcanoplot(
    reformat_de_table(dds_subset, x, annotation = features_db),
    x_cutoff = 0.25,
    main = paste(x[[2]], 'vs', x[[3]]))

  save_svg_plot(
    p,
    file.path(
      plasma_path, 'volcano-plot',
      paste0('de_volcano-plot_heparin-effect_additive-model_', x[[2]], '_vs_', x[[3]], '.svg')),
    width = 8, height = 8)
}))

## sars results models comparisons  ----

comparisons <- list(
  c("sars", "inactivated_sars", "no_sars"),
  c("sars", "infectious_sars", "no_sars"),
  c("sars", "infectious_sars", "inactivated_sars")
)

### heparin effect study, simple vs additive models, no thresholds ----
invisible(lapply(comparisons, compare_models))

### heparin effect study, simple vs additive models, PVAL < 0.05 ----
invisible(lapply(comparisons, compare_models, padj_threshold = 0.05))

# session info ----

rm(list=ls())

gc()

sessioninfo::session_info()

