# Load libraries ----

## GSEA analysis
# SOURCE: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)

## Data manipulation
library(data.table)

## Visualization
library(ggplot2)

## Set seed
set.seed(123)

# Objectives ----

# - Pairwise comparison of all SARS-CoV-2 categories (no_sars vs
#   inactivated_sars vs infectious_sars) using GSEA with MSigDB and GO
# - Decide if inactivated_sars and infectious_sars can be merged as a single
#   group.

# User-defined parameters ----

## Global variables ----
onlyhuman <- FALSE
genestag  <- if (onlyhuman) '_onlyhuman' else ''

# Adjusted p-value used to pre-filter the DE table (NULL = no filtering).
PVAL <- NULL # 0.1

# Genes used to build the ranking: 'ALL', 'UP-GENES' or 'DOWN-GENES'.
EXP <- 'ALL'
stopifnot(EXP %in% c('ALL', 'UP-GENES', 'DOWN-GENES'))

# GSEA parameters.
GSEA_PVAL   <- 0.05 # MSigDB collections
GSEGO_PVAL  <- 0.1  # GO sub-ontologies
MIN_GS_SIZE <- 10
MAX_GS_SIZE <- 800

# Tags used in output folder names (same naming as previous versions).
PVALtag <- if (is.null(PVAL)) '' else paste0('_PVAL', PVAL)
EXPtag  <- if (EXP == 'ALL') '' else paste0('_', EXP)

# fgsea score type: one-tailed tests when only up- or down-regulated genes
# are ranked.
score_type <- switch(
  EXP,
  'ALL'        = 'std',
  'UP-GENES'   = 'pos',
  'DOWN-GENES' = 'neg')

## Gene-set libraries ----

# MSigDB collections to analyse (name used in output folders = collection).
MSIGDB_COLLECTIONS <- c(
  'HALLMARK' = 'H') #,
  #'C2' = 'C2',
  #'C3' = 'C3',
  #'C5' = 'C5',
  #'C7' = 'C7')

## Paths ----

input_path <- file.path(
  project_path,
  'output',
  '05-dge',
  paste0('sars-effect-analysis', genestag),
  'de-tables')
stopifnot(dir.exists(input_path))

output_path <- file.path(
  project_path,
  'output',
  '06-gsea',
  paste0('sars-effect-analysis', genestag, PVALtag, EXPtag),
  'GSEA')
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## Comparisons ----

comparisons <- list(
  inactivated_sars_vs_no_sars         = c('sars', 'inactivated_sars', 'no_sars'),
  infectious_sars_vs_no_sars          = c('sars', 'infectious_sars', 'no_sars'),
  infectious_sars_vs_inactivated_sars = c('sars', 'infectious_sars', 'inactivated_sars'))

## Custom functions ----

source("r-utils/r-gsea.R")
source("r-utils/r-plots.R")
source("r-utils/r-tables.R")

# Main processing ----

## Load gene-set libraries ----
# Goal: load every gene-set library once, before looping over the p-value
# thresholds, since they do not depend on the DE tables.

msigdb_libs <- lapply(
  names(MSIGDB_COLLECTIONS),
  function(mdb) load_msigdb(MSIGDB_COLLECTIONS[[mdb]], prefix = mdb)) |>
  setNames(names(MSIGDB_COLLECTIONS))

dt_GO <- as.data.table(GO.db::GOTERM)[, .(go_id, Ontology, Definition)] |>
  unique()
setkey(dt_GO, 'go_id')

## Run GSEA per p-value threshold ----

### Load DE tables ----
# Goal: read every comparison once, apply the user-defined filters and
# build its ranked gene list.

de_tables <- setNames(
  lapply(names(comparisons), function(cmp) read_de_table(
    file.path(input_path, paste0('de_sars-effect_', cmp, '.csv')),
    padj_threshold = PVAL,
    direction = EXP,
    only_human = TRUE,
    drop_na_log2fc = TRUE,
    cmp_id = cmp)),
  names(comparisons))

# Skip comparisons left empty after filtering.
de_tables <- de_tables[vapply(de_tables, nrow, integer(1)) > 0]
if (length(de_tables) == 0) next

ranked_lists <- lapply(de_tables, make_ranked_list)
# Also evaluated with: 
#ranked_lists <- lapply(de_tables, make_ranked_list, rank_by = 'stat')

### GSEA w/ msigdbr libraries ----
# Goal: run fgsea-based GSEA on each MSigDB collection and, for every
# significant term, plot the log2FC of its leading-edge genes.

for (mdb in names(msigdb_libs)[1]) {
  
  msigdb       <- msigdb_libs[[mdb]]
  db_path      <- file.path(output_path, paste0('GSEA_msigdbr_w', mdb))
  heatmap_path <- file.path(db_path, 'individual_heatmaps')
  dir.create(heatmap_path, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(
    msigdb$genesets,
    file = file.path(output_path, paste0(mdb, '_definitions.tsv')),
    sep = '\t')
  
  for (cmp in names(ranked_lists)) {
    
    gse <- GSEA(
      ranked_lists[[cmp]],
      pvalueCutoff = GSEA_PVAL,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      pAdjustMethod = 'fdr',
      by = 'fgsea',
      scoreType = score_type,
      TERM2GENE = msigdb$term2gene,
      seed = TRUE)
    if (is_empty_result(gse)) next
    
    # Keep leading-edge genes as Ensembl IDs before converting them to SYMBOL.
    res_ens  <- as.data.frame(gse)
    core_ids <- setNames(res_ens$core_enrichment, res_ens$ID)
    rm(res_ens)
    
    # Change gene names from ENSEMBL to SYMBOL and add term descriptions.
    gse <- DOSE::setReadable(gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL')
    df_gse <- as.data.frame(gse)
    df_gse$Description <- msigdb$term2desc[df_gse$ID, description]
    
    write.csv(
      df_gse,
      file = file.path(db_path, paste0('table_sars-effect_', cmp, '.csv')),
      row.names = FALSE)
    
    plot_gsea_results(
      gse,
      ranked_lists[[cmp]],
      title = paste0(mdb, ': ', cmp),
      out_dir = db_path,
      cmp_id = cmp,
      analysis_tag = 'sars-effect',
      emap_tag = 'sars-terms_')
    
    # Heatmaps (one per term).
    for (term in names(core_ids)) {
      p <- do_heatmap_dge(
        de_tables[[cmp]],
        split_core_enrichment(core_ids[[term]]),
        subtitle = gsub('_', ' ', term),
        title = 'SARS effect analysis')
      
      file_base <- file.path(
        heatmap_path,
        paste0('heatmap_sars-terms_sars-effect_', cmp, '_', term))
      
      save_svg_plot(
        p,
        paste0(file_base, '.svg'),
        width = 4, height = 3 + nrow(p$data) * 0.3)
      
      fwrite(p$data, paste0(file_base, '.tsv'), sep = '\t')
    }
  }
}

### gseGO w/ GO.db terms ----
# Goal: run GSEA on each GO sub-ontology (BP, MF, CC) separately and add
# the GO term definitions to the output tables.

fwrite(
  dt_GO,
  file = file.path(output_path, 'GO_definitions.tsv'),
  sep = '\t')

for (cmp in names(ranked_lists)) {
  for (goset in c('BP', 'MF', 'CC')) {
    
    go_path <- file.path(output_path, paste0('gseGO_w', goset))
    dir.create(go_path, recursive = TRUE, showWarnings = FALSE)
    
    gse <- gseGO(
      geneList = ranked_lists[[cmp]],
      ont = goset,
      keyType = 'ENSEMBL',
      OrgDb = org.Hs.eg.db,
      pvalueCutoff = GSEGO_PVAL,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      pAdjustMethod = 'fdr',
      scoreType = score_type,
      verbose = TRUE,
      seed = TRUE)
    if (is_empty_result(gse)) next
    
    # Change gene names from ENSEMBL to SYMBOL and add GO definitions.
    gse <- DOSE::setReadable(gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL')
    df_gse <- as.data.frame(gse)
    df_gse$Definition <- dt_GO[df_gse$ID, Definition]
    
    write.csv(
      df_gse,
      file = file.path(go_path, paste0('table_sars-effect_', cmp, '.csv')),
      row.names = FALSE)
    
    plot_gsea_results(
      gse,
      ranked_lists[[cmp]],
      title = paste0('gseGO:', goset, ' | ', cmp),
      out_dir = go_path,
      cmp_id = cmp,
      analysis_tag = 'sars-effect',
      dotplot_height = 14)
  }
}

# Session info ----

rm(list = ls())
gc()
sessioninfo::session_info()
