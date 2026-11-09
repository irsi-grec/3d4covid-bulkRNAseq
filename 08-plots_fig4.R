# Load libraries ----

## Differential expression
library(DESeq2)
library(limma)

## Data manipulation
library(data.table)

## Visualization
library(ggplot2)
library(pheatmap)

## Set seed
set.seed(123)

## Custom functions
source("r-utils/r-io.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")
source("r-utils/r-overview.R")

# Objectives ----

# Figure 4: effect of plasma (plasma_covid vs plasma_healthy).
# - PCA of the 36 plasma samples and of the 12 plasma samples w/o SARS-CoV-2.
# - GSEA lollipop plots (HALLMARK and GO).
# - log2FC heatmaps of selected pathways.
# - Per-sample expression heatmaps of selected pathways (DESeq2 normalized
#   and vst data, the latter w/ the SARS-CoV-2 effect removed).
#
# Where two options exist, both are run and labelled:
# - [SUBMITTED] reproduces the version sent for review.
# - [REVISED]   corrected/alternative version.

# User-defined parameters ----

## Global variables ----

# GSEA run to plot (must match the folders written by 06-1).
PVAL <- NULL # 0.1
EXP  <- 'ALL' # 'UP-GENES' or 'DOWN-GENES'

PVALtag <- if (is.null(PVAL)) '' else paste0('_PVAL', PVAL)
EXPtag  <- if (EXP == 'ALL') '' else paste0('_', EXP)

MSIGDB_COLLECTION <- 'HALLMARK'
GO_SETS <- c('BP' = 8.5, 'MF' = 10, 'CC' = 8.5) # GO sub-ontology = plot width

N_TOP_PC1 <- 30  # top PC1 loading genes plotted in the heatmap
PADJ_CUT  <- 0.05

## Paths ----

de_path <- file.path(
  project_path,
  'output',
  '05-dge',
  'plasma-effect-analysis',
  'de-tables')

gsea_path <- file.path(
  project_path,
  'output',
  '06-gsea')

output_path <- file.path(
  project_path,
  'output',
  '08-plots',
  'fig4')
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

subset_path <- file.path(
  output_path,
  'subset_plasma-no_sars')
dir.create(subset_path, recursive = TRUE, showWarnings = FALSE)

## Comparisons ----

comparisons <- list(
  plasma_healthy_vs_no_plasma    = c('plasma', 'plasma_healthy', 'no_plasma'),
  plasma_covid_vs_no_plasma      = c('plasma', 'plasma_covid', 'no_plasma'),
  plasma_covid_vs_plasma_healthy = c('plasma', 'plasma_covid', 'plasma_healthy'))

TARGET_CMP <- 'plasma_covid_vs_plasma_healthy'

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

select_pathways <- list(
  `ecm` =
    c("COL16A1", "COL27A1", "ELN", "MMP2", "MMP14",
      "TIMP1", "TIMP3", "PLAU", "ADAM19", "ADAMTSL4"),
  `actin` =
    c("ACTR2", "ARF6", "ARHGAP1", "ARHGAP6", "ARHGAP10",
      "CDC42", "LMOD1", "MYH10", "RHOBTB1", "RHOBTB2",
      "SNTB2", "CDH4", "CLDN11", "CLDN23", "JAM3",
      "OCLN", "TJP2"),
  `adhesion` =
    c("ICAM3", "ICAM4", "IGFBP4", "IGFBP6", "PTPRU",
      "PTPN21", "PTPN22", "VCAM1"),
  `angiogen` =
    c("ADM", "CCND2", "CCNG1", "CTNNB1", "CCNL1",
      "FGF2", "FGFR3", "NOTCH3", "PDGFA", "PDGFB",
      "PDGFC", "PDGFD", "PDGFRA", "PDGFRB", "PDGFRL",
      "SLIT2", "SNTA1", "SMAD7", "SPRY4", "SPRED3",
      "TGFBI", "VEGFA", "VEGFB", "VEGFC"),
  `tnfa` =
    c("CXCL10", "DDX58", # RIGI in recent annotations (ENSG00000107201)
      "CSF2", "DUSP2", "GPR183", "LIF", "NR4A3",
      "HES1", "IL6", "TNFSF9", "SIK1", "DUSP5",
      "CCNL1", "EGR2", "NR4A2", "KLF10", "VEGFA",
      "NR4A1", "SLC2A3", "BHLHE40", "PTX3", "KLF2",
      "IL6ST", "JUNB", "SGK1", "HBEGF", "KLF6",
      "CDKN1A", "DNAJB4", "NFAT5"),
  `tnfa.2` =
    c("CXCL10", "DDX58",
      "CSF2", "DUSP2", "GPR183", "LIF", "NR4A3",
      "HES1", "IL6", "TNFSF9", "SIK1", "DUSP5",
      "CCNL1", "EGR2", "NR4A2", "KLF10", "VEGFA",
      "NR4A1", "SLC2A3", "BHLHE40", "PTX3"),
  `tnfa.3` =
    c("DDX58",
      "HES1", "TNFSF9", "DUSP5", "CCNL1", "KLF10",
      "VEGFA", "BHLHE40", "IL6ST", "HBEGF", "KLF6",
      "CDKN1A", "DNAJB4", "NFAT5")
)

# Plots per pathway: file tag, subtitle and plot heights (in inches) of the
# log2FC heatmap, the DESeq2-normalized heatmap and the vst heatmap.
pathway_plots <- data.table(
  set         = c('ecm', 'actin', 'adhesion', 'angiogen', 'tnfa'),
  tag         = c('ecm_remodelling', 'actin_regulation', 'cell_adhesion',
                  'angiongenesis_and_vasculogenesis', 'TNFa_signalling'),
  subtitle    = c('ECM remodelling', 'Actin regulation and cell junction', 'Cell adhesion',
                  'angiogenesis and vasculogenesis', 'TNFa signalling'),
  dge_height  = c(5, 5, 5, 8, 20),
  norm_height = c(3, 4, 3, 5, 7),
  vst_height  = c(3, 4, 3, 5, 12))

## Functions ----

# Save a plasma PCA as <file_name>.svg.
save_plasma_pca <- function(pca_res, md, title, out_path = output_path, file_name, shape = NULL) {
  p <- do_plotPCA_from_res(
    pca_res,
    md,
    mytitle = title,
    group.by = 'plasma',
    ellipse.by = 'plasma',
    manual_colours = plasma_palette,
    shape = shape,
    label = 'id')

  save_svg_plot(p, file.path(out_path, paste0(file_name, '.svg')), width = 8, height = 8)
}

# Main processing ----

## PCA w/ 36 plasma samples ----
# Goal: plasma_healthy and plasma_covid samples, all SARS-CoV-2 conditions.

md_36 <- metadata[plasma != 'no_plasma']

dds_36 <- make_dds_subset(
  dds,
  metadata,
  samples = md_36$id,
  vars = c('batch', 'sars', 'plasma'),
  design = ~ sars + plasma)
dds_36 <- DESeq(dds_36, test = 'Wald')

norm_mat_36 <- counts(dds_36, normalized = TRUE)

### PCAs, vst-transformed ----
# Goal: run PCA on variance-stabilized data, selecting the top 500 most
# variable genes on the same scale.

vsd_36 <- DESeq2::vst(
  dds_36,        # normalized counts cannot be used as input
  blind = FALSE, # use the design to estimate the dispersion trend
  fitType = 'parametric',
  nsub = 1000)

pca_res <- run_pca_top_var(assay(vsd_36))

save_plasma_pca(
  pca_res, md_36[, .(id, batch, sars, plasma)],
  title = 'PCA on vst data for plasma effect analysis',
  out_path = file.path(output_path, subset_path),
  file_name = 'pca-prcomp_plasma_36samples',
  shape = 'sars')

# Genes contributing most to PC1 (absolute loading).
loadings_PC1 <- pca_res$rotation[, 1]
top_PC1 <- data.table(
  GENEID = names(loadings_PC1),
  loading = loadings_PC1)[order(-abs(loading))]
top_PC1[, GENELAB := features_db[.(top_PC1$GENEID), GENELAB]]

fwrite(
  top_PC1[seq_len(N_TOP_PC1)],
  file.path(output_path, paste0('pca-prcomp_plasma_36samples_PC1-top', N_TOP_PC1, '.tsv')),
  sep = '\t')

top_pc1_genes <- top_PC1[seq_len(N_TOP_PC1), GENELAB]

rm(pca_res, loadings_PC1, vsd_36)

### PCAs, DESeq2-normalized counts [SUBMITTED] ----
# Goal: reproduce the submitted PCA (top 500 most variable genes on
# DESeq2-normalized counts, linear scale). Without a log transformation the
# variance is dominated by highly expressed genes.

save_plasma_pca(
  run_pca_top_var(norm_mat_36), md_36[, .(id, batch, plasma)],
  title = 'PCA on DESeq2 normalized data for plasma effect analysis',
  file_name = 'pca-prcomp_DESeqNorm_plasma_36samples')

### PCAs, log2-CPM values [REVISED] ----
# Goal: select the top 500 most variable genes and run the PCA on the same
# scale (TMM-normalized log2-CPM), so that gene selection is not driven by
# the mean-variance relationship of raw counts.

save_plasma_pca(
  run_pca_lcpm(counts(dds_36)), md_36[, .(id, batch, plasma)],
  title = 'PCA on log2-CPM data for plasma effect analysis',
  file_name = 'pca-prcomp_log2CPM_plasma_36samples')

## PCA w/ 12 plasma samples w/o SARS-CoV-2 ----

md_12 <- metadata[plasma != 'no_plasma' & sars == 'no_sars']

dds_12 <- make_dds_subset(
  dds,
  metadata,
  samples = md_12$id,
  vars = c('batch', 'plasma'),
  design = ~ plasma)
dds_12 <- DESeq(dds_12, test = 'Wald')

norm_mat_12 <- counts(dds_12, normalized = TRUE)

### PCAs, vst-transformed ----
# Note: the original call used the default ellipse.by = 'condition', a
# column not passed to the plot; ellipses are drawn by plasma.

vsd_12 <- DESeq2::vst(
  dds_12,
  blind = FALSE,
  fitType = 'parametric',
  nsub = 1000)

save_plasma_pca(
  run_pca_top_var(assay(vsd_12)), md_12[, .(id, batch, plasma)],
  title = 'PCA on vst data for plasma effect analysis',
  out_path = file.path(output_path, subset_path),
  file_name = 'pca-prcomp_plasma')

### PCAs, DESeq2-normalized counts [SUBMITTED] ----

save_plasma_pca(
  run_pca_top_var(norm_mat_12), md_12[, .(id, batch, plasma)],
  title = 'PCA on DESeq2 normalized data for plasma effect analysis',
  out_path = file.path(output_path, subset_path),
  file_name = 'pca-prcomp_DESeqNorm_plasma')

### PCAs, log2-CPM values [REVISED] ----

save_plasma_pca(
  run_pca_lcpm(counts(dds_12)), md_12[, .(id, batch, plasma)],
  title = 'PCA on log2-CPM data for plasma effect analysis',
  out_path = file.path(output_path, subset_path),
  file_name = 'pca-prcomp_log2CPM_plasma')

## GSEA lollipop plots ----
# Goal: plot the significant terms of plasma_covid vs plasma_healthy by
# GeneRatio (HALLMARK and every GO sub-ontology).

# GSEA result folders and plot widths (in inches).
gsea_runs   <- c(paste0('GSEA_msigdbr_w', MSIGDB_COLLECTION), paste0('gseGO_w', names(GO_SETS)))
gsea_widths <- c(7, unname(GO_SETS))

for (i in seq_along(gsea_runs)) {

  db_dir <- gsea_runs[[i]]

  df_gse <- read_gsea_table(
    gsea_path, 'plasma-effect', db_dir, TARGET_CMP,
    tag = paste0(PVALtag, EXPtag))
  if (is.null(df_gse) || nrow(df_gse) == 0) next

  # GO terms are labelled as '<Description>| <GO ID>'.
  if (startsWith(db_dir, 'gseGO')) df_gse[, ID := paste(Description, ID, sep = '| ')]
  df_gse <- reformat_data_dotplot(df_gse)
  
  df_gse_nes <- copy(df_gse)
  df_gse_nes[, ID := factor(as.character(ID), levels = as.character(ID)[order(-NES)])]
  
  p <- plot_gsea_lollipop(
    df_gse, x = 'GeneRatio',
    fixed_limits = FALSE)
  
  save_svg_plot(
    p, 
    file.path(output_path, paste0('plot_plasma-effect_', PVALtag, EXPtag, TARGET_CMP, '_', db_dir, '.svg')),
    width = gsea_widths[[i]], height = 5)
  
  p <- plot_gsea_lollipop(
    df_gse_nes, x = 'NES',
    x_max = 2.5, breaks_by = 0.5, label_x = 2,
    fixed_limits = FALSE)
  
  save_svg_plot(
    p, 
    file.path(output_path, paste0('plot_plasma-effect_', PVALtag, EXPtag, TARGET_CMP, '_', db_dir, '-NES.svg')),
    width = gsea_widths[[i]], height = 5)
}

## log2FC heatmaps of selected pathways ----

dge <- fread(file.path(de_path, paste0('de_plasma-effect_', TARGET_CMP, '.csv')))
dge[, cmp := factor(TARGET_CMP)]

for (i in seq_len(nrow(pathway_plots))) {
  pp <- pathway_plots[i]
  genes_set <- select_pathways[[pp$set]]

  p <- do_heatmap_dge_symbol(
    dge, genes_set, pp$subtitle,
    title = 'Plasma healthy vs plasma covid',
    padj_cut = PADJ_CUT)

  file_base <- file.path(
    output_path,
    paste0('heatmap_plasma-effect_', PVALtag, EXPtag, TARGET_CMP, '_dge-', pp$tag))

  save_svg_plot(p, paste0(file_base, '.svg'), width = 3.5, height = pp$dge_height)
  fwrite(dge[SYMBOL %in% genes_set], paste0(file_base, '.tsv'), sep = '\t')
}

## Per-sample expression heatmaps ----

### DESeq2 normalized data ----

for (i in seq_len(nrow(pathway_plots))) {
  pp <- pathway_plots[i]

  p <- do_heatmap_gex_per_sample(norm_mat_36, select_pathways[[pp$set]], md_36, features_db)

  save_pheatmap_svg(
    p,
    file.path(output_path, paste0('heatmap_plasma-effect_DESeqNorm_', PVALtag, EXPtag,
                                  TARGET_CMP, '-', pp$tag, '.svg')),
    width = 10, height = pp$norm_height)
}

### vst data w/o SARS-CoV-2 effect ----
# Goal: vst without design, then remove the SARS-CoV-2 effect with limma.
# - [SUBMITTED] removeBatchEffect(batch = sars), no design. Original file
#   names.
# - [REVISED] removeBatchEffect(batch = sars, design = ~ plasma), which
#   protects the plasma effect ('_design-plasma').

dds_36_null <- make_dds_subset(
  dds,
  metadata,
  samples = md_36$id,
  vars = c('batch', 'sars', 'plasma'),
  design = ~ 1)

vsd <- DESeq2::vst(
  dds_36_null,
  blind = FALSE,
  fitType = 'parametric',
  nsub = 1000)

expr_adj <- list(
  submitted = list(
    suffix = '',
    mat = limma::removeBatchEffect(
      assay(vsd),
      batch = vsd$sars)),
  revised = list(
    suffix = '_design-plasma',
    mat = limma::removeBatchEffect(
      assay(vsd),
      batch = vsd$sars,
      design = model.matrix(~ plasma, data = colData(vsd)))))

message('Max. absolute difference between both corrections: ',
        signif(max(abs(expr_adj$submitted$mat - expr_adj$revised$mat)), 3))

rm(dds_36_null, vsd)

# Reduced TNFa gene sets and top PC1 genes.
extra_heatmaps <- list(
  list(genes = select_pathways$tnfa.2, tag = 'TNFa_signalling-20genes', height = 5),
  list(genes = select_pathways$tnfa.3, tag = 'TNFa_signalling-14genes', height = 4),
  list(genes = top_pc1_genes,          tag = paste0('PCA-top', N_TOP_PC1), height = 6))

for (ea in expr_adj) {

  for (i in seq_len(nrow(pathway_plots))) {
    pp <- pathway_plots[i]

    p <- do_heatmap_gex_per_sample(ea$mat, select_pathways[[pp$set]], md_36, features_db)

    save_pheatmap_svg(
      p,
      file.path(output_path, paste0('heatmap_plasma-effect_', PVALtag, EXPtag,
                                    TARGET_CMP, '-', pp$tag, ea$suffix, '.svg')),
      width = 10, height = pp$vst_height)
  }

  for (eh in extra_heatmaps) {
    p <- do_heatmap_gex_per_sample(ea$mat, eh$genes, md_36, features_db)

    save_pheatmap_svg(
      p,
      file.path(output_path, paste0('heatmap_plasma-effect_', PVALtag, EXPtag,
                                    TARGET_CMP, '-', eh$tag, ea$suffix, '.svg')),
      width = 10, height = eh$height)
  }
}

## Heatmaps w/ 12 plasma samples w/o SARS-CoV-2 ----
# Goal: same pathway heatmaps for the subset plasma_healthy + plasma_covid
# without SARS-CoV-2 (design ~ plasma). No SARS-CoV-2 effect has to be removed
# because all samples are no_sars. Output folder: 'subset_plasma-no_sars'.

### log2FC heatmaps of selected pathways ----
# DE results of this subset (plasma_covid vs plasma_healthy, design ~ plasma).

dge_12 <- reformat_de_table(
  dds_12,
  comparisons[[TARGET_CMP]],
  annotation = features_db)
dge_12[, cmp := factor(TARGET_CMP)]

fwrite(
  dge_12,
  file.path(subset_path, paste0('de_plasma-effect_', TARGET_CMP, '_subset_plasma-no_sars.tsv')),
  sep = '\t')

for (i in seq_len(nrow(pathway_plots))) {
  pp <- pathway_plots[i]
  genes_set <- select_pathways[[pp$set]]
  
  p <- do_heatmap_dge_symbol(
    dge_12, genes_set, pp$subtitle,
    title = 'Plasma healthy vs plasma covid (no SARS-CoV-2)',
    padj_cut = PADJ_CUT)
  
  file_base <- file.path(
    subset_path,
    paste0('heatmap_plasma-effect_', PVALtag, EXPtag, TARGET_CMP, '_dge-', pp$tag))
  
  save_svg_plot(p, paste0(file_base, '.svg'), width = 3.5, height = pp$dge_height)
  fwrite(dge_12[SYMBOL %in% genes_set], paste0(file_base, '.tsv'), sep = '\t')
}

### Per-sample expression heatmaps ----
# DESeq2 normalized and vst data (vst w/ design ~ plasma, blind = FALSE).

expr_12 <- list(
  list(mat = counts(dds_12, normalized = TRUE), prefix = 'heatmap_plasma-effect_DESeqNorm_', height = 'norm_height'),
  list(mat = assay(vsd_12),                     prefix = 'heatmap_plasma-effect_',          height = 'vst_height'))

# Reduced TNFa gene sets (the top PC1 genes come from the 36-sample PCA).
extra_heatmaps_12 <- list(
  list(genes = select_pathways$tnfa.2, tag = 'TNFa_signalling-20genes', height = 5),
  list(genes = select_pathways$tnfa.3, tag = 'TNFa_signalling-14genes', height = 4))

for (ex in expr_12) {
  
  for (i in seq_len(nrow(pathway_plots))) {
    pp <- pathway_plots[i]
    
    p <- do_heatmap_gex_per_sample(ex$mat, select_pathways[[pp$set]], md_12, features_db)
    
    save_pheatmap_svg(
      p,
      file.path(subset_path, paste0(ex$prefix, PVALtag, EXPtag, TARGET_CMP, '-', pp$tag, '.svg')),
      width = 10, height = pp[[ex$height]])
  }
  
  for (eh in extra_heatmaps_12) {
    p <- do_heatmap_gex_per_sample(ex$mat, eh$genes, md_12, features_db)
    
    save_pheatmap_svg(
      p,
      file.path(subset_path, paste0(ex$prefix, PVALtag, EXPtag, TARGET_CMP, '-', eh$tag, '.svg')),
      width = 10, height = eh$height)
  }
}

rm(dds_12, vsd_12, md_12, dge_12, expr_12)

# Session info ----

rm(list = ls())
gc()
sessioninfo::session_info()
