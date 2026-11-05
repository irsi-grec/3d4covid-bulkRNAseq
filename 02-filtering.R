# Load libraries ----

## Differential expression
library(DESeq2)

## Data manipulation
library(data.table)
library(tidyr)

## Visualization
library(ggplot2)

## Set seed
set.seed(123)

# Objectives ----

# - Remove genes based on:
#   - Having 0s in all datasets.
#   - Minimum number of datasets with the given gene.

# User-defined parameters ----

## Paths ----

data_path <- file.path(
  project_path, 
  'output',
  '01-data-formating')

output_path <- file.path(
  project_path,
  'output',
  '02-filtering')

dir.create(output_path, recursive = TRUE)

## Custom functions ----

source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

## DGEList object

dge <- readRDS(
  file.path(data_path, paste0("DGEList_obj_", mapper, '.rds')))

## Metadata

metadata <- readRDS(
  file.path(data_path, paste0("metadata_", mapper,".rds")))

## Palettes

source("r-utils/r-palettes.R")

# Overview raw counts ----

# Perform log2CPM transformation.
lcpm <- cpm(
  dge,
  log = TRUE)

# Barplot total counts

# Add raw feature counts to metadata.
metadata[, nCounts.raw := colSums(dge$counts)[samples_id]]

# Set ascending order by condition.
metadata$id <- factor(
  metadata$id,
  levels = metadata$id[order(metadata$condition, decreasing = FALSE)])

# Generate barplot and save.
ggplot(metadata, aes(x = id, y = nCounts.raw, fill = condition)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  scale_fill_brewer(palette = 'Paired') +
  theme_minimal() +
  labs(x = '', y = 'Raw counts')

ggsave(
  file.path(output_path, 'barplot_raw_counts-by-condition.svg'),
  device = grDevices::svg,
  units = 'in',
  height = 4,
  width = 12,
  bg = 'white')

# Generate a table with the gene counts per sample and genes, and also
# genes information.
dt_counts <- data.table(
  GENEID = row.names(dge$counts),
  dge$counts
)

# Add gene information.
dt_counts <- merge(
  dt_counts,
  dge$genes,
  by = 'GENEID',
  all.x = TRUE)

fwrite(
  dt_counts,
  file.path(output_path, "dge_raw_counts.tsv"),
  sep = "\t")

remove(dt_counts)

## Gene variance ----

plot(rowMeans(dge$counts), rowSds(dge$counts), 
     main='Raw counts: mean vs sd')

plot(rowMeans(lcpm), rowSds(lcpm), 
     main='Raw logCPM: mean vs sd')

qplot(rowMeans(dge$counts), rowSds(dge$counts), log="xy") +
  geom_abline(intercept = 0, slope = 1, col="red") + #Poisson
  geom_smooth(col="orange") +
  theme_minimal()

## logCPM density before filtering ----

# Reformat data.

lcpm_long <- pivot_longer(
  as.data.frame(lcpm),
  cols = samples_id,
  names_to = 'samples',
  values_to = 'lcpm') |> 
  as.data.table()

lcpm_long$samples <- factor(
  lcpm_long$samples,
  levels = samples_id)

# Do density plot.

ggplot(
  lcpm_long, aes(x = lcpm, color = samples)) + 
  facet_grid(scales = "free_y") +
  geom_density(adjust=1.5) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x=element_blank()) + 
  labs(x = 'log2(CPM)') +
  ylim(0,1)

ggsave(
  file.path(output_path, 'density_plot_raw_counts.svg'),
  device = grDevices::svg,
  width = 6,
  height = 5,
  units = 'in',
  limitsize = FALSE,
  bg = 'white')

remove(lcpm_long)

# Filtering of genes ----

# Remove genes with 0 counts in all samples.
keep_genes <- rowSums(dge$counts) > 0
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
rm(keep_genes)

# TMM normalization.
dge <- calcNormFactors(
  dge,
  method = "TMM")

# Filtering of CPM using the experimental design.
design <- model.matrix(
  ~ condition,
  data = metadata)

keep_cpm <- filterByExpr(
  dge,
  design)

dge <- dge[keep_cpm, , keep.lib.sizes = FALSE]

# Save filtered matrix ----

# DGE object

dge_filter <- DGEList(
  counts = dge$counts,
  samples = metadata,
  group = metadata$condition,
  genes = dge$genes
)

rm(dge, keep_cpm)

saveRDS(
  dge_filter,
  file.path(output_path, paste0("DGEList_obj_", mapper,".filter.rds")))

write.csv(
  dge_filter$counts,
  file = file.path(output_path, "filter_matrix.csv"))

# DESeq object

dds <- DESeqDataSetFromMatrix(
  countData = dge_filter$counts,
  colData = metadata[, .(condition, batch, plasma, heparin, sars)],
  design = ~ 1)

saveRDS(
  dds,
  file.path(output_path, paste0("DESeq2_obj_", mapper, '.filter.rds')))

# Overview filtered counts ----

# Perform log2CPM transformation.
lcpm_filter <- cpm(
  dge_filter,
  log = TRUE) 

# Barplot total counts

# Add filtered feature counts to metadata.
metadata[, nCounts.filt := colSums(dge_filter$counts)[samples_id]]

# Set ascending order by condition.
metadata$id <- factor(
  metadata$id,
  levels = metadata$id[order(metadata$condition, decreasing = FALSE)])

# Generate barplot and save.
ggplot(metadata, aes(x = id, y = nCounts.filt, fill = condition)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  scale_fill_brewer(palette = 'Paired') +
  theme_minimal() +
  labs(x = '', y = 'Filter counts')

ggsave(
  file.path(output_path, 'barplot_filter_counts-by-condition.svg'), 
  device = grDevices::svg,
  units = 'in', height = 4, width = 12, bg = 'white')

# Barplot counts loss
metadata[, filtering.loss := round(((nCounts.raw - nCounts.filt) * 100) / nCounts.raw, 3)]

ggplot(metadata, aes(
  x = id, y = filtering.loss, fill = condition)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  scale_fill_brewer(palette = 'Paired') +
  theme_minimal() +
  labs(x = '', y = 'Percentage of filtered reads')

### Gene variance ----

# tentative plots
plot(rowMeans(dge_filter$counts), rowSds(dge_filter$counts), 
     main='Filter counts: mean vs sd')

plot(rowMeans(lcpm_filter), rowSds(lcpm_filter), 
     main='Filter logCPM: mean vs sd')

qplot(rowMeans(dge_filter$counts), rowSds(dge_filter$counts), log="xy") +
  geom_abline(intercept = 0, slope = 1, col="red") + #Poisson
  geom_smooth(col="orange") +
  theme_minimal()

## logCPM density after filtering ----

# Reformat data.

lcpm_filter_long <- pivot_longer(
  as.data.frame(lcpm_filter),
  cols = samples_id,
  names_to = 'samples',
  values_to = 'lcpm') |> 
  as.data.table()

lcpm_filter_long$samples <- factor(
  lcpm_filter_long$samples,
  levels = metadata$id)

# Do density plot.

ggplot(
  lcpm_filter_long, aes(x = lcpm, color = samples)) + 
  facet_grid(scales = "free_y") +  # facet_grid(rows = "samples", scales = "free_y") + 
  geom_density(adjust=1.5) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x=element_blank()) + 
  labs(x = 'log2(CPM)') +
  ylim(0,1)

ggsave(
  file.path(output_path, 'density_plot_filter_counts.svg'),
  device = grDevices::svg,
  width = 6, height = 5, units = 'in', limitsize = FALSE, bg = 'white')

rm(lcpm, lcpm_filter_long)

## PCAs ----

#### log2-CPM values - NEW ----
# Goal: select the top 500 most variable genes and run the PCA on the same
# scale (TMM-normalized log2-CPM), so that gene selection is not driven by
# the mean-variance relationship of raw counts.

pca_res <- run_pca_lcpm(dge_filter$counts)

save_pca_plots(
  pca_res, metadata,
  groups = c('condition', 'batch', 'plasma', 'heparin', 'sars'),
  titles = 'PCA on filtered data',
  out_dir = output_path,
  file_prefix = 'pca-prcomp_filter_cpm')

#### log2-CPM values - OLD ----

# Select top 500 most variable genes.
rv <- rowVars(dge_filter$counts)
select <- order(rv, decreasing = TRUE)[1:500]

pca_res <- prcomp(t(lcpm_filter[select,]))
remove(rv, select)

do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on filtered data',
  group.by = 'condition',
  manual_colours = conditions_palette,
  label = 'id')

ggsave(
  file.path(output_path, 'pca-prcomp_filter_counts_condition.svg'), 
  device = grDevices::svg,
  units = 'in', height = 8, width = 8, bg = 'white')

do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on filtered data',
  group.by = 'batch',
  ellipse.by = 'batch',
  manual_colours = c(`b1` = "#239B56", `b2` = "#3498DB", `b3` = "#E74C3C"),
  label = 'id')

ggsave(
  file.path(output_path, 'pca-prcomp_filter_counts_batch.svg'),
  device = grDevices::svg,
  units = 'in', height = 8, width = 8, bg = 'white')

do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on filtered data',
  group.by = 'plasma',
  ellipse.by = 'plasma',
  manual_colours = plasma_palette,
  label = 'id')

ggsave(
  file.path(output_path, 'pca-prcomp_filter_counts_plasma.svg'), 
  device = grDevices::svg,
  units = 'in', height = 8, width = 8, bg = 'white')

do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on filtered data',
  group.by = 'heparin',
  ellipse.by = 'heparin',
  manual_colours = heparin_palette,
  label = 'id')

ggsave(
  file.path(output_path, 'pca-prcomp_filter_counts_heparin.svg'), 
  device = grDevices::svg,
  units = 'in', height = 8, width = 8, bg = 'white')

do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on filtered data',
  group.by = 'sars',
  ellipse.by = 'sars',
  manual_colours = sars_palette,
  label = 'id')

ggsave(
  file.path(output_path, 'pca-prcomp_filter_counts_sars.svg'),
  device = grDevices::svg,
  units = 'in', height = 8, width = 8, bg = 'white')

rm(pca_res)

# Save filtered matrix metadata ----

saveRDS(
  metadata,
  file.path(output_path, paste0("metadata_", mapper,".filter.rds")))

write.csv(
  metadata,
  file = file.path(output_path, paste0("metadata_", mapper,".filter.csv")),
  row.names = FALSE)

# Session info ----
sessioninfo::session_info()

