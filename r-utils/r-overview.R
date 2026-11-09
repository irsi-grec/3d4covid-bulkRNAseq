# Shared helpers for the overview (07-*) and figure (08-*) scripts.
# Requires: data.table, ggplot2, DESeq2, pheatmap (loaded by the caller).

# DESeq2 subsets ----

# Build a DESeqDataSet for a subset of samples, with colData matched to the
# count matrix by sample id and unused factor levels removed.
make_dds_subset <- function(dds, metadata, samples, vars, design) {
  samples <- as.character(samples)
  md <- metadata[match(samples, as.character(id))]
  stopifnot(!anyNA(md$id))

  col_data <- data.frame(md[, vars, with = FALSE], row.names = samples)

  dds_sub <- DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(dds)[, samples],
    colData = col_data,
    design = design)

  for (v in vars) dds_sub[[v]] <- droplevels(dds_sub[[v]])
  dds_sub
}

# GSEA lollipop plots ----

# Read a GSEA table written by the 06-* scripts (NULL if it does not exist).
# `tag` is the PVALtag + EXPtag suffix of the analysis folder.
read_gsea_table <- function(gsea_path, effect, db_dir, cmp_id, tag = '') {
  gsea_file <- file.path(
    gsea_path,
    paste0(effect, '-analysis', tag),
    'GSEA',
    db_dir,
    paste0('table_', effect, '_', cmp_id, '.csv'))

  if (!file.exists(gsea_file)) {
    warning('GSEA table not found: ', gsea_file)
    return(NULL)
  }
  fread(gsea_file)
}

# Add Count (leading-edge genes), signed GeneRatio and Regulation, and order
# terms: activated first, then suppressed.
reformat_data_dotplot <- function(df_gsea) {
  dt <- as.data.table(df_gsea)
  dt[, Count := lengths(strsplit(core_enrichment, '/', fixed = TRUE))]
  dt[, Regulation := fifelse(enrichmentScore < 0, 'Suppressed', 'Activated')]
  dt[, GeneRatio := fifelse(Regulation == 'Suppressed', -Count / setSize, Count / setSize)]
  setorder(dt, Regulation, -GeneRatio)
  dt[, ID := factor(ID, levels = unique(ID))]
  dt
}

# Terms whose value falls outside the default x range of the lollipop plot.
gsea_out_of_range <- function(dt, x, x_max) {
  dt[abs(get(x)) > x_max, .(ID, value = get(x), p.adjust)]
}

# Lollipop plot of GSEA terms.
# fixed_limits = TRUE  -> x range fixed at +/- x_max (submitted version):
#                         terms outside the range are dropped by ggplot.
# fixed_limits = FALSE -> x range widened when needed, so no term is dropped.
plot_gsea_lollipop <- function(dt, x, x_max = 0.55, breaks_by = 0.25, label_x = 0.45, fixed_limits = TRUE) {
  x_lim  <- if (fixed_limits) x_max else max(x_max, max(abs(dt[[x]])) * 1.05)
  breaks <- seq(-floor(x_lim / breaks_by) * breaks_by, floor(x_lim / breaks_by) * breaks_by, breaks_by)
  y_top  <- nlevels(dt$ID) + 0.5

  ggplot(dt, aes(x = .data[[x]], y = ID, size = Count, color = p.adjust)) +
    geom_segment(aes(x = 0, xend = .data[[x]], y = ID, yend = ID),
                 color = 'gray60', linewidth = 0.4, linetype = 'solid') +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 'solid', color = 'black', linewidth = 0.8) +
    scale_color_gradient(low = '#3945A5', high = '#D1393A') +
    scale_x_continuous(limits = c(-x_lim, x_lim), breaks = breaks) +
    annotate('text', x = -label_x, y = y_top, label = 'Supressed', color = 'black', size = 4.5, hjust = 0) +
    annotate('text', x = label_x, y = y_top, label = 'Activated', color = 'black', size = 4.5, hjust = 1) +
    labs(title = '', x = x, y = '') +
    theme_minimal()
}

# log2FC heatmap of selected genes (by SYMBOL) ----

# One column per comparison (`cmp`); genes sorted by padj; '*' if padj < padj_cut.
# Genes missing from the DE table (e.g. alias symbols) are reported.
do_heatmap_dge_symbol <- function(de_table, target_genes, subtitle, title = NULL, padj_cut = 0.05) {

  missing <- setdiff(target_genes, de_table$SYMBOL)
  if (length(missing) > 0) {
    message(subtitle, ': genes not found in DE table: ', paste(missing, collapse = ', '))
  }

  dt <- de_table[SYMBOL %in% target_genes, .(SYMBOL, log2FoldChange, padj, cmp)]
  dt[is.na(padj), padj := max(de_table$padj, na.rm = TRUE)]

  ggplot(dt, aes(cmp, reorder(SYMBOL, padj, decreasing = TRUE), fill = log2FoldChange)) +
    geom_tile(color = 'white', lwd = 1.5, linetype = 1) +
    geom_text(aes(label = fifelse(padj < padj_cut, '*', '')), color = 'black', size = 5) +
    scale_fill_gradient2(low = '#3954A5', mid = '#FFFFFF', high = '#D1393A', midpoint = 0) +
    scale_y_discrete(drop = FALSE) +
    xlab('Comparison') +
    ylab('Genes') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(angle = 0, hjust = 0.5, size = 14),
      axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
    ggtitle(title, subtitle = subtitle)
}

# Per-sample expression heatmaps ----

# Map gene labels (SYMBOL, or GENEID when there is no symbol) to the gene IDs
# present in a matrix. Missing and duplicated labels are reported.
labels_to_ids <- function(gene_labels, available_ids, annotation) {
  map <- annotation[GENEID %in% available_ids & GENELAB %in% gene_labels, .(GENEID, GENELAB)]

  missing <- setdiff(gene_labels, map$GENELAB)
  if (length(missing) > 0) message('Genes not found in matrix: ', paste(missing, collapse = ', '))

  dup <- unique(map[duplicated(GENELAB), GENELAB])
  if (length(dup) > 0) message('Duplicated symbols, keeping the first ID: ', paste(dup, collapse = ', '))

  map[!duplicated(GENELAB)]
}

# Heatmap of z-scored expression (per gene) for selected genes and samples.
# Genes w/o variance (z-score NaN) are removed.
do_heatmap_gex_per_sample <- function(mat, target_genes, target_metadata, annotation) {
  samples <- as.character(target_metadata$id)
  map <- labels_to_ids(target_genes, rownames(mat), annotation)

  mat <- as.matrix(mat[map$GENEID, samples, drop = FALSE])
  rownames(mat) <- map$GENELAB

  mat <- t(scale(t(mat)))
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]

  annotation_col <- data.frame(
    plasma = target_metadata$plasma,
    sars = target_metadata$sars,
    row.names = samples)

  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    color = colorRampPalette(c('#3954A5', 'white', '#D1393A'))(50),
    annotation_col = annotation_col,
    silent = TRUE)
}

# Save a pheatmap object as SVG.
save_pheatmap_svg <- function(p, file, width = 8, height = 8) {
  svg(file, width = width, height = height, bg = 'white')
  on.exit(dev.off())
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  invisible(file)
}
