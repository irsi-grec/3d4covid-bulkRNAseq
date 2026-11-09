# Define custom theme ----

theme_myplot <- function() {
  theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    )
}

# Generate boxplots ----

do_boxplot <- function(dt, x = "sample", y = "count", x_fill = x,
                       mytitle = NULL, xtitle = NULL, ytitle = NULL) {
  
  stopifnot(is.data.table(dt))
  stopifnot(all(c(x, y, x_fill) %in% names(dt)))
  
  y_data <- dt[[y]]
  y_data[is.infinite(y_data)] <- NA
  
  ggplot(dt, aes(x = .data[[x]], y = y_data, fill = .data[[x_fill]])) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Paired") +
    theme_myplot() +
    labs(title = mytitle, x = xtitle, y = ytitle)
}

# PCA ----

# Auxilary funtion for PCA

plot_PCA_generic <- function(pca_dt, group.by = "condition", label = NULL,
                             shape = NULL, manual_colours = NULL,
                             mytitle = NULL, ellipse.by = group.by) {
  
  stopifnot(is.data.table(pca_dt))
  stopifnot(all(c("PC1", "PC2") %in% names(pca_dt)))
  stopifnot(group.by %in% names(pca_dt))
  
  if (is.null(mytitle)) {
    mytitle <- paste("PCA colored by", group.by)
  }
  
  # Point aesthetics: colour by group.by and, optionally, shape by `shape`.
  point_aes <- if (is.null(shape)) {
    aes(color = .data[[group.by]])
  } else {
    aes(color = .data[[group.by]], shape = .data[[shape]])
  }
  
  p <- ggplot(pca_dt, aes(x = PC1, y = PC2)) +
    geom_point(point_aes, size = 2) +
    stat_ellipse(aes(color = .data[[ellipse.by]]), type = "norm") +
    theme_myplot() +
    labs(
      x = sprintf("PC1: %.2f%% variance", pca_dt$var.explained[1] * 100),
      y = sprintf("PC2: %.2f%% variance", pca_dt$var.explained[2] * 100),
      title = mytitle,
      color = group.by
    )
  
  if (!is.null(manual_colours)) {
    p <- p + scale_color_manual(values = manual_colours)
  }
  
  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = .data[[label]], color = .data[[group.by]]),
      size = 3, box.padding = 0
    )
  }
  
  p
}

# Main funtion
do_plotPCA_from_res <- function(res, metadata, group.by = "condition",
                                ellipse.by = "condition", label = NULL,
                                shape = NULL, manual_colours = NULL,
                                mytitle = NULL) {
  
  pca_dt <- data.table(
    PC1 = res$x[,1],
    PC2 = res$x[,2],
    var.explained = summary(res)$importance[2,],
    metadata
  )
  
  plot_PCA_generic(
    pca_dt,
    group.by = group.by,
    label = label,
    shape = shape,
    manual_colours = manual_colours,
    mytitle = mytitle,
    ellipse.by = ellipse.by
  )
}

# PCA on the top `ntop` most variable genes of a (genes x samples) matrix.
# Gene selection and PCA are done on the same scale.
run_pca_top_var <- function(mat, ntop = 500) {
  rv <- matrixStats::rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  prcomp(t(mat[select, ]))
}

# PCA on TMM-normalized log2-CPM values (edgeR), so that gene selection is not
# driven by the mean-variance relationship of raw counts.
run_pca_lcpm <- function(counts, ntop = 500, prior.count = 2) {
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = 'TMM')
  lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
  run_pca_top_var(lcpm, ntop)
}

# PCA on variance-stabilized data (DESeq2::vst, blind = FALSE).
run_pca_vst <- function(dds, ntop = 500) {
  run_pca_top_var(SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE)), ntop)
}

# Save one PCA plot per grouping variable as <out_dir>/<file_prefix>_<group>.svg.
# `titles`: a single title for all plots, or a vector named by group.
save_pca_plots <- function(pca_res, metadata, groups, titles, out_dir, file_prefix,
                           palettes = pca_palettes, label = 'id',
                           width = 8, height = 8) {
  for (g in groups) {
    p <- do_plotPCA_from_res(
      pca_res,
      metadata,
      mytitle = if (length(titles) > 1) titles[[g]] else titles,
      group.by = g,
      ellipse.by = g,
      manual_colours = palettes[[g]],
      label = label)
    
    save_svg_plot(p, file.path(out_dir, paste0(file_prefix, '_', g, '.svg')), width, height)
  }
  invisible(NULL)
}

# Autoplot ----

do_autoplot <- function(pca = pca.res, data = metadata, colour = 'condition', label = NULL, shape = NULL, manual_colours = NULL, title = NULL){
  myp <- autoplot(pca,
                  data = data,
                  colour = colour,
                  shape = shape,
                  size = 3) + 
    scale_fill_brewer(palette = 'Paired')  + 
    theme_minimal() #theme_ipsum()
  
  if(!is.null(manual_colours)){
    myp <- myp + scale_color_manual(values = manual_colours)
  }
  
  if(!is.null(label)){
    myp <- myp + geom_text_repel(aes(x=PC1, y=PC2, label=get(label)), box.padding = 0, size = rel(3))
  }
  
  if(!is.null(title)){
    myp <- myp + ggtitle(title)
  }
  
  return(myp)
}

# Volcano plot ----


do_volcanoplot <- function(dt, x = "log2FoldChange", y = "padj",
                           Label = "SYMBOL", top = 10, topnames = NULL,
                           x_cutoff = log2(1.5), y_cutoff = 0.05,
                           main = NULL, xlab = "Log2 Fold Change",
                           ylab = "-Log10(P)", control_guides = "NO_SITE",
                           label_fontsize = 4, ...) {
  
  stopifnot(is.data.table(dt))
  stopifnot(all(c(Label, x, y) %in% names(dt)))
  
  gg <- dt[, .(Label = get(Label), x = get(x), y = get(y))]
  
  gg[, hit := fifelse(grepl(control_guides, Label), "ctrl",
                      fifelse(x > x_cutoff & y < y_cutoff, "up",
                              fifelse(x < -x_cutoff & y < y_cutoff, "down", "a")))]
  
  gg[, y_plot := -log10(fifelse(y == 0, min(y[y > 0], na.rm = TRUE), y))]
  
  # ranking
  gg[, rank := frank(y, ties.method = "random")]
  gg[, id_label := fifelse(Label %in% topnames | rank <= top, Label, "")]
  
  colors <- c(a = "gray80", ctrl = "#D4A017", up = "#e41a1c", down = "#377eb8")
  
  p <- ggplot(gg, aes(x = x, y = y_plot, color = hit)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted") +
    geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted") +
    scale_color_manual(values = colors) +
    theme_myplot() +
    labs(title = main, x = xlab, y = ylab)
  
  if (any(gg$id_label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = gg[id_label != ""],
      aes(label = id_label),
      size = label_fontsize
    )
  }
  
  p
}


# plotDists ----

#The distance plot calculates the Euclidean distance between the expression values 
#for each individual sample. This is a way to identify which samples are most 
#closely related in terms of their gene expression patterns.

plot_distance <- function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  p <- pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
  p
}

# Variable genes heatmap ----

variable_gene_heatmap <- function(vsd.obj, num_genes = 500, annotation,
                                  title = "") {
  
  stabilized <- assay(vsd.obj)
  vars <- matrixStats::rowVars(stabilized)
  
  top_genes <- stabilized[order(vars, decreasing = TRUE)[1:num_genes], ]
  top_genes <- sweep(top_genes, 1, rowMeans(top_genes), "-")
  
  setkey(annotation, "GENEID")
  gene_names <- annotation[rownames(top_genes), GENELAB]
  rownames(top_genes) <- gene_names
  
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  
  lim <- max(abs(top_genes))
  n_col <- 256
  
  pheatmap::pheatmap(
    top_genes,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n_col),
    breaks = seq(-lim, lim, length.out = n_col + 1),
    annotation_col = coldata,
    fontsize_col = 8,
    fontsize_row = 250 / num_genes,
    border_color = NA,
    main = title
  )
}

# Heatmap of log2FC for the leading-edge genes of one term ----

# Genes are matched by Ensembl ID, so no manual SYMBOL aliasing is needed
# between EnsDb.Hsapiens.v86 and org.Hs.eg.db.
do_heatmap_dge <- function(de_table, target_ids, subtitle, title = NULL) {
  
  dt <- de_table[.(target_ids), .(id, SYMBOL, log2FoldChange, padj, cmp), nomatch = NULL]
  dt[, gene_label := fifelse(is.na(SYMBOL) | SYMBOL == '', id, SYMBOL)]
  dt[, gene_label := make.unique(gene_label)]
  dt[, gene_label := factor(gene_label, levels = sort(unique(gene_label), decreasing = TRUE))]
  
  ggplot(dt, aes(x = cmp, y = gene_label, fill = log2FoldChange)) +
    geom_tile(color = 'white', lwd = 1.5, linetype = 1) +
    geom_text(aes(label = fifelse(!is.na(padj) & padj < 0.05, '*', '')),
              color = 'black', size = 5) +
    scale_fill_gradient2(low = '#3954A5', mid = '#FFFFFF', high = '#D1393A', midpoint = 0) +
    scale_y_discrete(drop = FALSE) +
    labs(x = 'Comparison', y = 'Genes') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(angle = 0, hjust = 0.5, size = 14),
      axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    ggtitle(title, subtitle = subtitle)
}

# Save plots as SVG ----

# Save a plot as SVG. Works with ggplot objects and with base/grid graphics
# expressions (e.g. gplots::venn(x)), which are evaluated once the device is open.
save_svg_plot <- function(p, file, width = 8, height = 8) {
  svg(file, width = width, height = height, bg = 'white')
  on.exit(dev.off())
  if (inherits(p, 'ggplot')) print(p) else p
  invisible(file)
}