# Load libraries ----

## Differential expression
library(DESeq2)

## Data manipulation
library(data.table)

## Visualization
library(ggplot2)
library(ggrepel)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")
source("r-utils/r-overview.R")

# Objectives ----

# Figure 2: effect of SARS-CoV-2 (no_plasma samples).
# - Sample distances, variable genes heatmap and PCA.
# - log2FC scatter plot of inactivated vs infectious SARS-CoV-2 (gene
#   categories and log2FC correlation).
# - GSEA (HALLMARK) lollipop plots by GeneRatio and NES.
# - log2FC heatmaps of selected pathways.
#
# Where two options exist, both are run and labelled:
# - [SUBMITTED] reproduces the version sent for review.
# - [REVISED]   corrected/alternative version.

# User-defined parameters ----

## Global variables ----

# GSEA run to plot (must match the folders written by 06-2).
PVAL <- NULL # 0.1
EXP  <- 'ALL' # 'UP-GENES' or 'DOWN-GENES'

PVALtag <- if (is.null(PVAL)) '' else paste0('_PVAL', PVAL)
EXPtag  <- if (EXP == 'ALL') '' else paste0('_', EXP)

MSIGDB_COLLECTION <- 'HALLMARK'

LFC_CUT  <- 0.25
PADJ_CUT <- 0.05

sars_palette2 <- c(
  `no_sars` = 'black',
  `inactivated_sars` = '#3498DB',
  `infectious_sars` = '#239B56')

## Paths ----

de_path <- file.path(
  project_path,
  'output',
  '05-dge',
  'sars-effect-analysis',
  'de-tables')

gsea_path <- file.path(
  project_path,
  'output',
  '06-gsea')

output_path <- file.path(
  project_path,
  'output',
  '08-plots',
  'fig2')
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## Comparisons ----

comparisons <- list(
  inactivated_sars_vs_no_sars         = c('sars', 'inactivated_sars', 'no_sars'),
  infectious_sars_vs_no_sars          = c('sars', 'infectious_sars', 'no_sars'),
  infectious_sars_vs_inactivated_sars = c('sars', 'infectious_sars', 'inactivated_sars'))

## Load data ----

inputs <- load_dge_inputs(project_path, mapper)
dds <- inputs$dds
metadata <- inputs$metadata
features_db <- inputs$features_db
setkey(features_db, 'GENEID')
rm(inputs)

## Load palettes ----

source("r-utils/r-palettes.R")

## Gene sets ----

# One log2FC heatmap per pathway: genes, comparisons shown and plot size.
heatmap_sets <- list(
  list(
    tag = 'dge-ecm', subtitle = 'ECM',
    cmps = names(comparisons)[1:2], width = 3.8, height = 6,
    genes = c("ADAM11", "ADAMTS17", "COL15A1", "COL1A1", "COL21A1",
              "COL4A3BP", "COL5A1", "ELN", "LAMA1", "LAMA4",
              "MMP15", "MMP19", "MMP23B", "MMP3")),
  list(
    tag = 'dge-junction', subtitle = 'Junction',
    cmps = names(comparisons)[1:2], width = 4, height = 8,
    genes = c("ACTN4", "CFL2", "FBN1", "FBN2", "LIMK1",
              "MYH10", "PFN4", "RAC3", "RHOF", "TPM1",
              "TPM2", "CDC42", "CDH13", "CLDN11", "CLDN12",
              "CLDN23", "ROCK2")),
  list(
    tag = 'dge-cell_adhesion', subtitle = 'Cell adhesion',
    cmps = names(comparisons)[1:2], width = 4, height = 4,
    genes = c("ANGPTL1", "ITGA1", "ITGA9", "ITGB1", "ITGB1BP2",
              "PODN", "VCAM1")),
  list(
    tag = 'dge-angiongenesis_and_vasculogenesis', subtitle = 'Ang & Vasc',
    cmps = names(comparisons)[1:2], width = 4, height = 12,
    genes = c("ADM", "APOE", "CXCL10", "EGF", "FGF1",
              "FGF14", "FGFR4", "FLNA", "FST", "FSTL1",
              "HIF3A", "IGFBP7", "IL6", "PDGFA", "PDGFC",
              "PDGFD", "PDGFRA", "PDGFRB", "PLAU", "PTK2B",
              "RHOC", "SERPINE1", "SLIT3", "SMAD6", "SMAD9",
              "THPO", "VEGFA", "VEGFB", "VEGFC", "VEGFD",
              "TGFBI", "TGFBR1")),
  list(
    tag = 'dge-TNFa_signalling', subtitle = 'TNFa signalling',
    cmps = names(comparisons), width = 4, height = 12,
    genes = c("AREG", "ATF3", "BCL3", "BHLHE40", "BIRC3",
              "CCL20", "CCNL1", "CD69", "CD83", "CDKN1A",
              "CEBPB", "CSF2", "CXCL3", "DUSP1", "DUSP2",
              "DUSP5", "EGR2", "EGR3", "ETS2", "F2RL1",
              "F3", "FJX1", "FOSB", "FOSL1", "FOSL2",
              "FUT4", "GPR183", "HES1", "ID2", "IL18",
              "IL1B", "IL23A", "IL6", "IRF1", "IRS2",
              "JUNB", "KLF10", "KLF2", "KLF4", "KLF9",
              "LIF", "MAFF", "MAP3K8", "MCL1", "MXD1",
              "NFIL3", "NFKBIA", "NR4A1", "NR4A2", "NR4A3",
              "PDE4B", "PER1", "PHLDA2", "PLAU", "PLAUR",
              "PLEK", "PTGS2", "RELB", "SERPINB2", "SERPINB8",
              "SERPINE1", "SGK1", "SIK1", "SLC2A3", "SPHK1",
              "SQSTM1", "TIPARP", "TNFAIP3", "TNFSF9", "TUBB2A",
              "VEGFA", "ZC3H12A", "ZFP36")),
  list(
    tag = 'dge-IFNa_signalling', subtitle = 'IFNa signalling',
    cmps = names(comparisons)[1:2], width = 3.5, height = 9,
    genes = c("C1S", "CASP8", "CCRL2", "CXCL10", "DDX60",
              "DHX58", "EPSTI1", "HERC6", "IFI35", "IFI44",
              "IFIT2", "IFITM1", "IFITM3", "IL7", "IRF1",
              "ISG20", "MOV10", "OAS1", "OASL", "PARP9",
              "PARP12", "PSMB9", "RSAD2", "RTP4", "SAMD9L",
              "SELL", "TXNIP",
              "WARS")) # synonym WARS1
)

## Functions ----

# Save the SARS-CoV-2 PCA with both palettes: <prefix>_sars.svg and
# <prefix>_sars_v2.svg.
save_sars_pca <- function(pca_res, title, prefix) {
  palettes <- list(sars_palette, sars_palette2)
  suffixes <- c('', '_v2')

  for (i in seq_along(palettes)) {
    p <- do_plotPCA_from_res(
      pca_res,
      metadata_subset[, .(id, batch, sars)],
      mytitle = title,
      group.by = 'sars',
      ellipse.by = 'sars',
      manual_colours = palettes[[i]],
      label = 'id')

    save_svg_plot(p, file.path(output_path, paste0(prefix, '_sars', suffixes[i], '.svg')), width = 8, height = 8)
  }
}

# Gene categories of the scatter plot.
# [SUBMITTED] 'lfc': categories defined only by |log2FC| >= LFC_CUT in each
#   comparison (genes already filtered by padj < PADJ_CUT in any comparison).
classify_lfc <- function(dt) {
  dt[, fcase(
    padj.inac_vs_no_sars > PADJ_CUT & padj.inf_vs_no_sars > PADJ_CUT,
      'non-significant',
    abs(log2FoldChange.inac_vs_no_sars) < LFC_CUT & abs(log2FoldChange.inf_vs_no_sars) < LFC_CUT,
      'non-significant',
    (log2FoldChange.inac_vs_no_sars >= LFC_CUT & log2FoldChange.inf_vs_no_sars >= LFC_CUT) |
      (log2FoldChange.inac_vs_no_sars <= -LFC_CUT & log2FoldChange.inf_vs_no_sars <= -LFC_CUT),
      'shared-gene',
    abs(log2FoldChange.inf_vs_no_sars) < LFC_CUT,
    'inactivated sars specific',
    abs(log2FoldChange.inac_vs_no_sars) < LFC_CUT,
    'infectious sars specific',
    default = 'divergent genes')]
}

# [REVISED] 'sig': a gene changes in a comparison if padj < PADJ_CUT and
#   |log2FC| >= LFC_CUT. 'specific' also requires a significant difference in
#   the direct comparison infectious vs inactivated; genes changing in one
#   comparison only, without that difference, are 'ambiguous'.
classify_sig <- function(dt) {
  sig_inac <- dt$padj.inac_vs_no_sars < PADJ_CUT & abs(dt$log2FoldChange.inac_vs_no_sars) >= LFC_CUT
  sig_inf  <- dt$padj.inf_vs_no_sars  < PADJ_CUT & abs(dt$log2FoldChange.inf_vs_no_sars)  >= LFC_CUT
  #sig_diff <- dt$padj.inf_vs_inac < PADJ_CUT
  same_dir <- sign(dt$log2FoldChange.inac_vs_no_sars) == sign(dt$log2FoldChange.inf_vs_no_sars)

  fcase(
    sig_inac & sig_inf & same_dir,  'shared-gene',
    sig_inac & sig_inf,             'divergent genes',
    #sig_inac & !sig_inf & sig_diff, 'inactivated sars specific',
    #!sig_inac & sig_inf & sig_diff, 'infectious sars specific',
    sig_inac & !sig_inf,            'inactivated sars specific',
    !sig_inac & sig_inf,            'infectious sars specific',
    sig_inac | sig_inf,             'ambiguous',
    default =                       'non-significant')
}

# Pearson correlation (r, r², 95% CI, p-value) between both log2FC, plus
# Spearman's rho as a robustness check against outliers.
cor_lfc <- function(dt, gene_set) {
  x <- dt$log2FoldChange.inac_vs_no_sars
  y <- dt$log2FoldChange.inf_vs_no_sars
  ct <- cor.test(x, y, method = 'pearson')

  data.table(
    gene_set     = gene_set,
    n            = length(x),
    r            = unname(ct$estimate),
    r2           = unname(ct$estimate)^2,
    ci95_low     = ct$conf.int[1],
    ci95_high    = ct$conf.int[2],
    p_value      = ct$p.value,
    spearman_rho = cor(x, y, method = 'spearman'))
}

# Main processing ----

## Subset samples ----
# Goal: keep only no_plasma samples and run DESeq2 (design ~ sars).

metadata_subset <- metadata[plasma == 'no_plasma']

dds <- make_dds_subset(
  dds,
  metadata,
  samples = metadata_subset$id,
  vars = c('batch', 'sars'),
  design = ~ sars)

dds <- DESeq(dds, test = 'Wald')
resultsNames(dds)

## Overview w/ vst data ----

vsd <- DESeq2::vst(
  dds,           # normalized counts cannot be used as input
  blind = FALSE, # use the design to estimate the dispersion trend
  fitType = 'parametric',
  nsub = 1000)

### distance plot ----

vsd_labeled <- vsd
colnames(vsd_labeled) <- paste(metadata_subset$condition, metadata_subset$id, sep = '|')

png(file = file.path(output_path, 'distance_plot.png'),
    units = 'in', height = 8, width = 8, bg = 'white', res = 1100)
plot_distance(vsd_labeled)
dev.off()

### variable genes heatmap ----

png(file = file.path(output_path, 'heatmap.png'),
    units = 'in', height = 8, width = 8, bg = 'white', res = 1100)
variable_gene_heatmap(
  vsd_labeled,
  num_genes = 50,
  annotation = features_db)
dev.off()

rm(vsd_labeled)

### PCAs, vst-transformed ----
# Goal: run PCA on variance-stabilized data, selecting the top 500 most
# variable genes on the same scale.

save_sars_pca(
  run_pca_top_var(assay(vsd)),
  title = 'PCA on vst data for sars effect analysis',
  prefix = 'pca-prcomp')

### PCAs, DESeq2-normalized counts [SUBMITTED] ----
# Goal: reproduce the submitted PCA (top 500 most variable genes on
# DESeq2-normalized counts, linear scale). Without a log transformation the
# variance is dominated by highly expressed genes.

save_sars_pca(
  run_pca_top_var(counts(dds, normalized = TRUE)),
  title = 'PCA on DESeq2 normalized data for sars effect analysis',
  prefix = 'pca-prcomp_DESeqNorm')

### PCAs, log2-CPM values [REVISED] ----
# Goal: select the top 500 most variable genes and run the PCA on the same
# scale (TMM-normalized log2-CPM), so that gene selection is not driven by
# the mean-variance relationship of raw counts.

save_sars_pca(
  run_pca_lcpm(counts(dds)),
  title = 'PCA on log2-CPM data for sars effect analysis',
  prefix = 'pca-prcomp_log2CPM')

## Load DE tables ----

de_tables <- rbindlist(lapply(names(comparisons), function(cmp_id) {
  dt <- fread(file.path(de_path, paste0('de_sars-effect_', cmp_id, '.csv')))
  dt[, cmp := cmp_id]
  dt
}))
de_tables[, cmp := factor(cmp, levels = names(comparisons))]

## Scatter plot inactivated vs infectious ----
# Goal: compare the log2FC of inactivated and infectious SARS-CoV-2 vs
# no_sars for the genes significant in any of both comparisons.

### Missing padj overview ----
# padj is NA when DESeq2 removes a gene by independent filtering (low mean
# counts; pvalue available) or flags it as a Cook's outlier (pvalue also NA).
# These genes are not considered significant.

scatter_cmps <- names(comparisons)[1:2]

na_summary <- de_tables[cmp %in% scatter_cmps, .(
  n_genes        = .N,
  padj_NA        = sum(is.na(padj)),
  pvalue_NA      = sum(is.na(pvalue)),
  indep_filtered = sum(!is.na(pvalue) & is.na(padj)),
  significant    = sum(padj < PADJ_CUT, na.rm = TRUE)),
  by = cmp]
print(na_summary)

### Build comparison table ----

diff_cmp <- dcast(
  de_tables[cmp %in% scatter_cmps][, cmp := as.character(cmp)],
  id + SYMBOL ~ cmp,
  value.var = c('padj', 'log2FoldChange'))

setnames(
  diff_cmp,
  c('padj_inactivated_sars_vs_no_sars', 'padj_infectious_sars_vs_no_sars',
    'log2FoldChange_inactivated_sars_vs_no_sars', 'log2FoldChange_infectious_sars_vs_no_sars'),
  c('padj.inac_vs_no_sars', 'padj.inf_vs_no_sars',
    'log2FoldChange.inac_vs_no_sars', 'log2FoldChange.inf_vs_no_sars'))

# padj of the direct comparison (used by the [REVISED] classification).
diff_cmp <- merge(
  diff_cmp,
  de_tables[cmp == 'infectious_sars_vs_inactivated_sars', .(id, padj.inf_vs_inac = padj)],
  by = 'id',
  all.x = TRUE)

# All genes w/ log2FC in both comparisons (reference set for the correlation).
diff_all <- na.omit(diff_cmp, cols = c('log2FoldChange.inac_vs_no_sars', 'log2FoldChange.inf_vs_no_sars'))

# Missing log2FC -> 0 (no change); missing padj -> 1 (not significant).
setnafill(diff_cmp, fill = 0, cols = c('log2FoldChange.inac_vs_no_sars', 'log2FoldChange.inf_vs_no_sars'))
setnafill(diff_cmp, fill = 1, cols = c('padj.inac_vs_no_sars', 'padj.inf_vs_no_sars', 'padj.inf_vs_inac'))

# Keep genes significant in any comparison.
diff_cmp <- diff_cmp[padj.inac_vs_no_sars < PADJ_CUT | padj.inf_vs_no_sars < PADJ_CUT]

# Label genes w/ large log2FC (Ensembl ID when there is no symbol).
diff_cmp[, GENELAB := fcoalesce(features_db$GENELAB[match(id, features_db$GENEID)], id)]
diff_cmp[, gene_label := fifelse(
  abs(log2FoldChange.inac_vs_no_sars) >= 2 | abs(log2FoldChange.inf_vs_no_sars) >= 2, GENELAB, '')]

### Correlation of log2FC ----

cor_summary <- rbind(
  cor_lfc(diff_cmp, 'significant in any comparison'), # scatter plot
  cor_lfc(diff_all, 'all genes w/ log2FC in both comparisons'))
print(cor_summary)

fwrite(
  cor_summary,
  file.path(output_path, 'scatterplot_sars_effect-inac.inf_vs_no-sars_correlation.tsv'),
  sep = '\t')

fwrite(
  na_summary,
  file.path(output_path, 'scatterplot_sars_effect-inac.inf_vs_no-sars_padj-NA.tsv'),
  sep = '\t')

### Gene categories [SUBMITTED] and [REVISED] ----
# Goal: plot and count the gene categories with both classifications. The
# [SUBMITTED] version keeps the original file names.

col_plot <- c(
  'non-significant' = 'grey',
  'shared-gene' = 'firebrick',
  'divergent genes' = 'black',
  'inactivated sars specific' = 'green4',
  'infectious sars specific' = 'dodgerblue',
  'ambiguous' = 'grey60')

classifications <- list(
  submitted = list(fun = classify_lfc, suffix = ''),
  revised   = list(fun = classify_sig, suffix = '_sig-based'))

for (cl in names(classifications)) {

  dt <- copy(diff_cmp)
  dt[, analysis := classifications[[cl]]$fun(dt)]

  scatter_summary <- dt[, .N, by = analysis][order(-N)]
  scatter_summary <- rbind(
    scatter_summary,
    data.table(analysis = 'total significant in any comparison', N = nrow(dt)))
  message('Gene categories [', toupper(cl), ']')
  print(scatter_summary)

  p <- ggplot(dt, aes(
    x = log2FoldChange.inac_vs_no_sars,
    y = log2FoldChange.inf_vs_no_sars,
    label = gene_label,
    col = analysis)) +
    geom_point(size = 0.5) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30) +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = 'dotted') +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 'dotted') +
    scale_color_manual(values = col_plot) +
    theme_minimal() +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 3))

  file_base <- file.path(
    output_path,
    paste0('scatterplot_sars_effect-inac.inf_vs_no-sars', classifications[[cl]]$suffix))

  save_svg_plot(p, paste0(file_base, '.svg'), width = 7, height = 7)
  fwrite(dt, paste0(file_base, '.tsv'), sep = '\t')
  fwrite(scatter_summary, paste0(file_base, '_summary.tsv'), sep = '\t')
}

## GSEA lollipop plots ----
# Goal: plot the significant HALLMARK terms of every comparison, by
# GeneRatio and by NES, with both x-axis options:
# - [SUBMITTED] fixed limits (+/-0.55 GeneRatio, +/-2.5 NES): terms outside
#   are dropped by ggplot. Original file names.
# - [REVISED] limits widened when needed ('_auto-limits').
# Terms outside the fixed limits are saved in '_out-of-range.tsv'.

for (cmp in names(comparisons)) {

  df_gse <- read_gsea_table(
    gsea_path,                                    
    'sars-effect',                                
    paste0('GSEA_msigdbr_w', MSIGDB_COLLECTION),  
    cmp,                                          
    tag = paste0(PVALtag, EXPtag))
  if (is.null(df_gse) || nrow(df_gse) == 0) next
  df_gse <- reformat_data_dotplot(df_gse)

  df_gse_nes <- copy(df_gse)
  df_gse_nes[, ID := factor(as.character(ID), levels = as.character(ID)[order(-NES)])]
  
  p <- plot_gsea_lollipop(
    df_gse, x = 'GeneRatio',
    fixed_limits = FALSE)
  save_svg_plot(
    p,
    file.path(output_path, paste0('plot_sars-effect_', cmp, '.svg')),
    width = 7, height = 5)
  
  p <- plot_gsea_lollipop(
    df_gse_nes, x = 'NES',
    fixed_limits = FALSE)
  save_svg_plot(
    p,
    file.path(output_path, paste0('plot_sars-effect_', cmp, '-NES', '.svg')),
    width = 7, height = 5)

}

## log2FC heatmaps of selected pathways ----

for (hs in heatmap_sets) {

  dge <- de_tables[cmp %in% hs$cmps]
  dge[, cmp := factor(as.character(cmp), levels = hs$cmps)]

  p <- do_heatmap_dge_symbol(
    dge, hs$genes, hs$subtitle,
    title = 'SARS effect analysis',
    padj_cut = PADJ_CUT)

  file_base <- file.path(
    output_path,
    paste0('heatmap_sars-effect_', PVALtag, EXPtag, hs$tag))

  save_svg_plot(p, paste0(file_base, '.svg'), width = hs$width, height = hs$height)
  fwrite(dge[SYMBOL %in% hs$genes], paste0(file_base, '.tsv'), sep = '\t')
}

# Session info ----

rm(list = ls())
gc()
sessioninfo::session_info()
