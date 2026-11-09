# Load the batch-corrected DESeq2 object, sample metadata and gene annotation
# (restricted to the genes present in the count matrix).
load_dge_inputs <- function(project_path, mapper = 'STAR') {
  data_path  <- file.path(project_path, 'output', '01-data-formating')
  input_path <- file.path(project_path, 'output', '03-batch-correction')

  dds <- readRDS(
    file.path(input_path, paste0('DESeq2_obj_', mapper, '.filter.combat.rds')))
  metadata <- readRDS(
    file.path(data_path, paste0('metadata_', mapper, '.rds')))
  features_db <- readRDS(file.path(data_path, 'genesDb.rds'))
  features_db <- features_db[row.names(SummarizedExperiment::assay(dds))]

  list(dds = dds, metadata = metadata, features_db = features_db)
}

# Count matrix of a DESeq2 object, optionally restricted to human genes.
# `genestag` is appended to output folder names.
select_genes <- function(dds, features_db, onlyhuman = FALSE) {
  mat <- SummarizedExperiment::assay(dds)
  if (onlyhuman) {
    mat <- mat[row.names(mat) %in% features_db[species == 'homo sapiens', GENEID], ]
  }
  list(mat = mat, genestag = if (onlyhuman) '_onlyhuman' else '')
}

# Create a base output folder and its sub-folders; returns the base path.
create_output_dirs <- function(base, subdirs = character()) {
  for (d in c(base, file.path(base, subdirs))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  invisible(base)
}

