# Build the ranked gene list required by GSEA: rank named by Ensembl ID,
# sorted in decreasing order.
#make_ranked_list <- function(de_table) {
#  stopifnot(anyDuplicated(de_table$id) == 0)
#  setNames(de_table$log2FoldChange, de_table$id) |>
#    sort(decreasing = TRUE)
#}
make_ranked_list <- function(de_table, rank_by = 'log2FoldChange') {
  stopifnot(anyDuplicated(de_table$id) == 0)
  de_table <- de_table[!is.na(get(rank_by))]
  setNames(de_table[[rank_by]], de_table$id) |>
    sort(decreasing = TRUE)
}

# Retrieve one MSigDB collection as TERM2GENE and TERM2DESC tables.
# msigdbr >= 10.0.0 renamed the `category` argument to `collection`.
load_msigdb <- function(collection, prefix) {
  gs <- if (packageVersion('msigdbr') >= '10.0.0') {
    msigdbr(species = 'Homo sapiens', collection = collection)
  } else {
    msigdbr(species = 'Homo sapiens', category = collection)
  }
  gs <- as.data.table(gs)
  gs[, gene_set := sub(paste0('^', prefix, '_'), '', gs_name)]

  term2gene <- unique(
    gs[!is.na(ensembl_gene) & ensembl_gene != '', .(gene_set, Geneid = ensembl_gene)])

  term2desc <- unique(
    gs[, .(gene_set, description = gs_description)],
    by = 'gene_set')
  setkey(term2desc, 'gene_set')

  list(genesets = gs, term2gene = term2gene, term2desc = term2desc)
}

# Split a clusterProfiler 'core_enrichment' string into a gene vector.
split_core_enrichment <- function(x) {
  strsplit(x, '/', fixed = TRUE) |>
    unlist() |>
    unique() |>
    sort()
}

# TRUE if a clusterProfiler result is NULL or has no significant terms.
is_empty_result <- function(gse) {
  is.null(gse) || nrow(as.data.frame(gse)) == 0
}

# GSEA plots from clusterProfiler ----

# Standard set of clusterProfiler plots for one GSEA result: gene-concept
# network, dotplot split by activated/suppressed, heatplot and enrichment map.
# Output files: <plot>_<emap_tag><analysis_tag>_<cmp_id>.svg
plot_gsea_results <- function(gse, ranked_list, title, out_dir, cmp_id,
                              analysis_tag,
                              emap_tag = '',
                              dotplot_height = 10,
                              emap_size = 10) {
  
  # Gene-concept network.
  # Disabled: gene-concept network not used.
  #p <- cnetplot(gse, categorySize = 'pvalue', color.params = list(foldChange = ranked_list)) +
  #  ggtitle(title)
  #save_svg_plot(
  #  p,
  #  file.path(out_dir, paste0('network_', analysis_tag, '_', cmp_id, '.svg')),
  #  width = 10, height = 10)
  
  # Dotplot, split by activated/suppressed pathways.
  p <- dotplot(gse, showCategory = 30, split = '.sign') +
    facet_grid(. ~ .sign) +
    ggtitle(title)
  save_svg_plot(
    p,
    file.path(out_dir, paste0('dotplot_', analysis_tag, '_', cmp_id, '.svg')),
    width = 10, height = dotplot_height)
  
  # Heatplot of the top terms.
  p <- heatplot(gse, foldChange = ranked_list, showCategory = 10) +
    ggtitle(title)
  save_svg_plot(
    p,
    file.path(out_dir, paste0('heatplot_', analysis_tag, '_', cmp_id, '.svg')),
    width = 25, height = 5)
  
  # Enrichment map (needs at least two terms).
  if (nrow(as.data.frame(gse)) >= 2) {
    p <- emapplot(pairwise_termsim(gse)) +
      ggtitle(title)
    save_svg_plot(
      p,
      file.path(out_dir, paste0('enrichmentmap_', emap_tag, analysis_tag, '_', cmp_id, '.svg')),
      width = emap_size, height = emap_size)
  }
}
