matrix2longerdt <- function(data = matrix, target_cols = cols_mat, factor_name = 'sample', value_name = 'count'){
  pivot_longer(
    data,
    cols = names(data[, ..target_cols]),
    names_to = factor_name,
    values_to = value_name) %>% 
    as.data.table()
}

reformat_de_table <- function(myde, mycmp, LFCCUT = NULL, PCUT = NULL, annotation = NULL){
  res <- results(myde, contrast = mycmp)
  
  dt <- as.data.table(res)
  dt[, id := row.names(res)]
  dt <- dt[, .(id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)]
  
  dt <- merge(dt, annotation, by.x = 'id', by.y = 'GENEID', all.x = TRUE)
  
  if(!is.null(LFCCUT)){
    dt <- dt[abs(log2FoldChange) > LFCCUT,]
  }
  
  if(!is.null(PCUT)){
    dt <- dt[padj < PCUT,]
  }
  setkey(dt, 'padj')
  dt
}

# Read DE tables ----

# Read a DE table written by reformat_de_table() and, optionally, filter rows,
# select columns and add a suffix to the value columns (all except id and
# SYMBOL), e.g. to merge several comparisons side by side.
read_de_table <- function(file,
                          padj_threshold = NULL,
                          log2fc_threshold = NULL,
                          direction = c('ALL', 'UP-GENES', 'DOWN-GENES'),
                          only_human = FALSE,
                          drop_na_log2fc = FALSE,
                          cols = NULL,
                          suffix = NULL,
                          cmp_id = NULL) {
  direction <- match.arg(direction)
  dt <- fread(file)
  
  if (!is.null(padj_threshold))   dt <- dt[padj < padj_threshold]
  if (!is.null(log2fc_threshold)) dt <- dt[abs(log2FoldChange) >= log2fc_threshold]
  if (direction == 'UP-GENES')    dt <- dt[log2FoldChange > 0]
  if (direction == 'DOWN-GENES')  dt <- dt[log2FoldChange < 0]
  if (only_human)                 dt <- dt[species == 'homo sapiens']
  if (drop_na_log2fc)             dt <- dt[!is.na(log2FoldChange)]
  
  if (!is.null(cols)) dt <- dt[, ..cols]
  
  if (!is.null(suffix)) {
    value_cols <- setdiff(names(dt), c('id', 'SYMBOL'))
    setnames(dt, value_cols, paste0(value_cols, '.', suffix))
  }
  
  if (!is.null(cmp_id)) dt[, cmp := cmp_id]
  
  setkey(dt, 'id')
  dt
}

# Build an output-file suffix from the thresholds used, e.g. '_PVAL0.05_LFC1'.
make_suffix <- function(padj_threshold = NULL, log2fc_threshold = NULL) {
  paste0(
    if (!is.null(padj_threshold))   paste0('_PVAL', padj_threshold),
    if (!is.null(log2fc_threshold)) paste0('_LFC',  log2fc_threshold),
    collapse = '')
}

