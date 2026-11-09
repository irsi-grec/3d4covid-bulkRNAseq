# Load libraries ----

## Differential expression
library(DESeq2)
library(gprofiler2)
library(edgeR)

## Data manipulation
library(data.table)
library(tidyr)

## Visualization
library(ggplot2)
library(RColorBrewer)
library(ggfortify)
library(ggrepel)
library(plotly)
library(htmlwidgets)

## Set seed
set.seed(123)

# Objectives ----

# - Check if there is batch effect in our samples.
# - Correct batch effect.

# User-defined parameters ----

## Global variables ----
onlyhuman <- FALSE

## Paths ----

data_path <- file.path(
  project_path, 
  'output',
  '01-data-formating')

input_path <- file.path(
  project_path,
  'output',
  '02-filtering')

output_path <- file.path(
  project_path,
  'output',
  '03-batch-correction')
dir.create(output_path, recursive = TRUE)

## Genes ----
features_db <- readRDS(file.path(data_path, "genesDb.rds"))

## Custom functions

source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

## DESeq2 object

dds <- readRDS(
  file.path(input_path, paste0("DESeq2_obj_", mapper, '.filter.rds')))

## Metadata

metadata <- readRDS(
  file.path(data_path, paste0("metadata_", mapper,".rds")))

## Palettes

source("r-utils/r-palettes.R")

# main ----

# Create a DGEList object.
# if required, select only human genes

if(onlyhuman){
  mat <- assay(dds)[row.names(assay(dds)) %in% features_db[species == 'homo sapiens', GENEID],]
  genestag <- '_onlyhuman'
  output_path <- file.path(paste0(output_path, genestag))
  dir.create(file.path(output_path, genestag), recursive = TRUE, showWarnings = FALSE)
} else {
  mat <- assay(dds)
  genestag <- ''
}

dge <- DGEList(
  counts = mat,
  samples = metadata,
  group = metadata$condition,
  genes = row.names(assay(dds))
)

# remove batch effect ----
design <- model.matrix(
  ~ plasma + heparin + sars,
  data = metadata)

dge <- edgeR::estimateGLMCommonDisp(
  dge,
  design = design)

adj <- sva::ComBat_seq(
  dge$counts,
  batch = metadata$batch,
  group = NULL, #metadata$condition,
  covar_mod = design,
  full_mod = TRUE)

# notice that ComBat_seq() return raw counts that should be treated as raw 
# counts # SOURCE: https://www.biostars.org/p/9573484/
identical(adj, dge$counts)
dge$counts <- adj

## save DGEList object ----
saveRDS(
  dge,
  file = file.path(
    output_path,
    paste0("DGEList_obj_", mapper,".filter.combat.rds")))

## save DESeq2 object ----

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = dge$counts,
  colData = metadata[, .(condition, batch, plasma, heparin, sars)],
  design = design)

rm(dge); gc()

mcols(dds) <- DataFrame(data.frame(gene=rownames(assay(dds))))

# Set factor levels and orders.

# set reference of levels
dds$batch <- factor(
  dds$batch,
  levels = c(
    "b1",
    "b2",
    "b3"))

dds$condition <- factor(
  dds$condition,
  levels = c(
    "control_mock",
    "heparin",
    "inactivated_control",
    "inactivated_healthy",
    "inactivated_heparin",
    "inactivated_covid19",
    "infected_control",
    "infected_healthy",
    "infected_heparin",
    "infected_covid19",
    "plasma_healthy",
    "plasma_covid19"))

dds$plasma <- factor(
  dds$plasma,
  levels = c(
    "no_plasma",
    "plasma_healthy",
    "plasma_covid"))

dds$heparin <- factor(
  dds$heparin,
  levels = c(
    "no_heparin",
    "heparin"))

dds$sars <- factor(
  dds$sars,
  levels = c(
    "no_sars",
    "inactivated_sars",
    "infectious_sars"))

saveRDS(
  dds,
  file = file.path(
    output_path,
    paste0("DESeq2_obj_", mapper,".filter.combat.rds")))

# plots ----

## gene variance ----

# tentative plots 
plot(rowMeans(assay(dds)), rowSds(assay(dds)), 
     main='Counts w/o batch: mean vs sd')

lcpm.dds <- cpm(assay(dds), log = TRUE)
plot(rowMeans(lcpm.dds), rowSds(lcpm.dds), 
     main='logCPM w/o batch: mean vs sd')

qplot(rowMeans(assay(dds)), rowSds(assay(dds)), log="xy") +
  geom_abline(intercept = 0, slope = 1, col="red") + #Poisson
  geom_smooth(col="orange") +
  theme_minimal()

rm(lcpm.dds)

# reported plot
log2mat <- log2(assay(dds) + 2)
df <- data.frame(
  Geneid = row.names(assay(dds)),
  SYMBOL = features_db[row.names(log2mat), .(SYMBOL)],
  rowMeans = rowMeans(log2mat),
  rowSds = rowSds(log2mat),
  rowVar = rowSds(log2mat)^2
)
df$label <- ifelse(df$rowVar >= 1.5, df$SYMBOL, '')

p <- ggplot(df, aes(rowMeans, rowSds)) +
  geom_point(size = 2) + 
  theme_minimal() +
  labs(x = "Mean of log2(expression + 2)",
       y = "Variance of log2(expression + 2)", # to avoid -Inf values when count = 0
       title = "Gene variance w/o batch") +
  ggrepel::geom_text_repel(aes(label = label), color = "black")
p
ggsave(file.path(output_path, 'scatterplot_gene-variance_woBatch.svg'), plot = p, units = 'in', height = 8, width = 8, bg = 'white')

rm(df)

## distance plot ----
png(file = file.path(output_path, 'distance_plot.png'),
    units = 'in', height = 8, width = 12, bg = 'white', res = 1100)
plot_distance(dds)
dev.off()

## variable genes heatmap ----
vsd <- DESeq2::vst(dds, blind = TRUE)
png(file = file.path(output_path, 'heatmap_vsd_blindTRUE.png'),
    units = 'in', height = 8, width = 12, bg = 'white', res = 1100)
variable_gene_heatmap(
  vsd,
  num_genes = 50,
  annotation = features_db)
dev.off()

vsd <- DESeq2::vst(dds, blind = FALSE)
png(file = file.path(output_path, 'heatmap_vsd_blindFALSE.png'),
    units = 'in', height = 8, width = 12, bg = 'white', res = 1100)
variable_gene_heatmap(
  vsd,
  num_genes = 50,
  annotation = features_db)
dev.off()

## pca ----

#### log2-CPM values - REVIEWED ----

pca_res <- run_pca_lcpm(assay(dds))

save_pca_plots(
  pca_res, metadata,
  groups = c('condition', 'heparin', 'sars', 'plasma', 'batch'),
  titles = c(
    condition = 'PCA on batch corrected data',
    heparin   = 'PCA on batch corrected data for heparin effect analysis',
    sars      = 'PCA on batch corrected data for sars effect analysis',
    plasma    = 'PCA on batch corrected data for plasma effect analysis',
    batch     = 'PCA on batch corrected data for batch effect analysis'),
  out_dir = output_path,
  file_prefix = 'pca-prcomp_cpm')

#### log2-CPM values - SUBMITTED ----

colnames(assay(dds))
rv <- rowVars(assay(dds))
select <- order(rv, decreasing = TRUE)[1:500]

pca_res <- prcomp(t(cpm(assay(dds), log = TRUE)[select,]))

p <- do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on batch corrected data',
  group.by = 'condition',
  manual_colours = conditions_palette,
  label = 'id')
p

svg(file.path(output_path, 'pca-prcomp_condition.svg'), width = 8, height = 8, bg = 'white')
print(p)
dev.off()

p <- do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on batch corrected data for heparin effect analysis',
  group.by = 'heparin',
  ellipse.by = 'heparin',
  manual_colours = heparin_palette,
  label = 'id')
p

svg(file.path(output_path, 'pca-prcomp_heparin.svg'), width = 8, height = 8, bg = 'white')
print(p)
dev.off()

p <- do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on batch corrected data for sars effect analysis',
  group.by = 'sars',
  ellipse.by = 'sars',
  manual_colours = sars_palette,
  label = 'id')
p

svg(file.path(output_path, 'pca-prcomp_sars.svg'), width = 8, height = 8, bg = 'white')
print(p)
dev.off()

p <- do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on batch corrected data for plasma effect analysis',
  group.by = 'plasma',
  ellipse.by = 'plasma',
  manual_colours = plasma_palette,
  label = 'id')
p

svg(file.path(output_path, 'pca-prcomp_plasma.svg'), width = 8, height = 8, bg = 'white')
print(p)
dev.off()

p <- do_plotPCA_from_res(
  pca_res,
  metadata,
  mytitle = 'PCA on batch corrected data for batch effect analysis',
  group.by = 'batch',
  ellipse.by = 'batch',
  manual_colours = c(`b1` = '#239B56', `b2` = '#3498DB', `b3` = '#E74C3C'),
  label = 'id')
p

svg(file.path(output_path, 'pca-prcomp_batch.svg'), width = 8, height = 8, bg = 'white')
print(p)
dev.off()

# session info ----
rm(list=ls())
gc()
sessioninfo::session_info()
