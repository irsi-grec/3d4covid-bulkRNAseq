# Load libraries ----

## Differential expression
library(edgeR)
library(limma)
library(DESeq2)

## Annotation
library(EnsDb.Hsapiens.v86)
library(gprofiler2)

## Data manipulation
library(data.table)
library(tidyr)

## Visualization
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(RColorBrewer)

## Set seed
set.seed(123)

# Objectives ----

# - Convert STAR count matrix into DGEList object
# - Add SARS-CoV-2 annotation into Human's gene annotation

# User-defined parameters ----

## Paths ----

# Output alignment with STAR is stored at star_path as indicated:
star_path <- file.path(
  project_path,
  'output',
  '00-alignment',
  'GRCh38.p14_MN908947.3',
  mapper
)

# Create output path
output_path <- file.path(
  project_path,
  'output',
  '01-data-formating'
)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

## Load count matrices ----
# Read the per-sample STAR gene-count files into a single DGEList,
# then split it into a plain counts matrix and a metadata table that will
# later be recombined with the experimental design information.

files <- list.files(
  path = star_path,
  pattern = '.*\\ReadsPerGene.out.tab',
  full.names = T)

# Extract the gene ID (1st column) and the correct counts (4th column).
dge_raw <- readDGE(
  files,
  columns = c(1, 4))

# STAR reports three summary rows before the gene counts.
# Remove them to keep only gene-level counts.
dge <- dge_raw[-seq_len(3), ]

# Extract sample IDs from STAR output file names.
# Example:
#   1_CM_ReadsPerGene.out.tab -> 1_CM
samples_id_prefix <- dge$samples$files |> 
  basename() |>
  gsub('^([^_]+_[^_]+).*', '\\1', x = _)

# Add metadata into the matrix
dge$samples$id <- samples_id_prefix
row.names(dge$samples) <- samples_id_prefix
colnames(dge$counts) <- samples_id_prefix

# Extract count matrix and add the genes symbols and set the order of columns.
features_matrix <- dge$counts |>
  as.data.table()

features_matrix <- features_matrix[, ..samples_id] |>
  as.data.frame()

rownames(features_matrix) <- row.names(dge$counts)

# Set up a data.table with meta data information from alignment.
features_metadata <- dge$samples |>
  as.data.table()

setkey(features_metadata, 'id')
features_metadata <- features_metadata[samples_id]

# remove unnecessary objects.
rm(star_path, samples_id_prefix, files, dge_raw, dge)

## Load metadata ----
# Goal: load the experimental design metadata, harmonize category labels,
# merge in the library-size information extracted from the alignment above,
# and encode all grouping variables as factors with explicit level order
# (required for correct contrasts in downstream differential expression).
metadata <- fread(
  file.path(project_path, 'data', 'samples_id_groups_batch.csv'))
setnames(metadata, 'group', 'condition')

# Re-format columns.
metadata$condition <- gsub(' ', '_', metadata$condition)
metadata[, conditions := paste0(
  metadata$plasma, ':', metadata$heparin, ':', metadata$sars)]

# Merge samples metadata and features information by samples' id.
metadata <- merge(
  metadata,
  features_metadata[, .(id, lib.size)],
  by = 'id')

# Set the levels' orders for every factor.
metadata$id <- factor(
  metadata$id,
  levels = samples_id)

metadata$batch <- factor(
  metadata$batch,
  levels = c('b1', 'b2', 'b3'))

metadata$sars <- factor(
  metadata$sars,
  levels = c('no_sars', 'inactivated_sars', 'infectious_sars'))

metadata$heparin <- factor(
  metadata$heparin,
  levels = c('no_heparin', 'heparin'))

metadata$plasma <- factor(
  metadata$plasma,
  levels = c('no_plasma', 'plasma_healthy', 'plasma_covid'))

metadata$condition <- factor(
  metadata$condition,
  levels = c(
    'control_mock',
    'inactivated_control',
    'infected_control',
    'heparin',
    'inactivated_heparin',
    'infected_heparin',
    'plasma_healthy',
    'inactivated_healthy',
    'infected_healthy',
    'plasma_covid19',
    'inactivated_covid19',
    'infected_covid19'))

# Re-order metadata's rows.
setkey(metadata, 'id')
metadata <- metadata[samples_id]

# set conditions metadata in order
metadata[, conditions := factor(
  paste0(metadata$plasma, ':', metadata$heparin, ':', metadata$sars),
  levels = unique(metadata$conditions))]

saveRDS(
  metadata,
  file = file.path(project_path, "data", "samples_id_groups_batch.rds")
)

# Main processing ----

## Create DGEList object ----
# Goal: assemble the final DGEList combining raw gene counts with the full
# sample metadata (experimental design + alignment QC info) built above.

dge <- DGEList(
  counts = features_matrix,
  samples = as.data.frame(metadata),
  group = metadata$condition
)

rm(features_metadata, features_matrix)

## Organize gene annotations ----
# Goal: retrieve human gene annotations (name, symbol, biotype) for every
# gene in the count matrix from Ensembl, and check for missing or
# duplicated annotation entries before attaching them to the DGEList.

# Extract genes annotation from Ensembldb
edb <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  keys = rownames(dge$counts),
  keytype = 'GENEID',
  columns = c('GENENAME', 'SYMBOL', 'GENEBIOTYPE')) |>
  as.data.table()

# Notice that the selected genes from Ensembldb may be fewer than all of our
# genes within the matrix.
stopifnot(
  nrow(edb) == nrow(dge$counts)
)

# Check duplicated genes' IDs (expected number = 0)
edb[duplicated(edb$GENEID), GENEID] |>
  sort() |> unique() |> length()

# Check duplicated genes' name (expected number != 0)
edb[duplicated(edb$GENENAME), GENENAME] |>
  sort() |> unique() |> length()

# One solution may be keeping only genes that share GeneId and Chr position
# (yet, many genes are removed with this option).
# Keep only genes that are in the matrix
# We used chromosome information, because duplicated genes can come from
# different regions

features_db <- merge(
  data.frame(GENEID = rownames(dge$counts)),
  edb,
  by = 'GENEID',
  all.x = TRUE) |>
  as.data.table()
setkey(features_db, 'GENEID')

# Check the number of genes w/o information.
features_db$GENENAME |>
  is.na() |>
  sum()

rm(edb)

## Add sars annotion ----
# Goal: the alignment reference combined the human genome with the
# SARS-CoV-2 genome (MN908947.3), so viral genes have no entry in
# EnsDb.Hsapiens.v86. This block parses the viral GTF annotation and merges
# it into the gene annotation table built above.

# Load SARS-CoV-2 annotation from GTF file.
sars_file <- fread(file.path(project_path, 'data', 'MN908947.3.gtf'))

# Add column names.
setnames(
  sars_file,
  names(sars_file),
  c('seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'))

# Keep only CDs.
sars_file <- sars_file[feature == 'CDS']

# Split attribute to extract gene name and id.
sars_file[, `:=`(
  GENENAME = sub('.*gene_name "([^"]+)".*', '\\1', attribute),
  SYMBOL   = sub('.*gene_name "([^"]+)".*', '\\1', attribute),  
  GENEID   = sub('.*gene_id "([^"]+)".*',   '\\1', attribute)
)]

# Merge Human and SARS-CoV-2 annotations.
features_db[, species := 'homo sapiens']
sars_file[, species := 'sars-cov-2']

features_db[sars_file, on = 'GENEID',
            `:=`(GENENAME = i.GENENAME,
                 SYMBOL   = i.SYMBOL,
                 species  = i.species)]

# Check, again, the number of genes w/o information.
features_db$GENENAME |>
  is.na() |>
  sum()

rm(sars_file)

## Add annotation into dge ----
# Goal: attach the completed gene annotation table (human + viral genes) to
# the DGEList, and add a display-friendly gene label column for plotting.

# Check concordant GENEIDs.
stopifnot(
  identical(features_db$GENEID, rownames(dge$counts))
)

# Sort genes table as dgeList object.
setkey(features_db, 'GENEID')
features_db <- features_db[row.names(dge$counts)]

# Add genes to DGEList.
table(rownames(dge$counts) == features_db$GENEID) # check the order
dge$genes <- features_db
names(dge)

# Add GENELAB column containing the gene name to our gene database
features_db[, GENELAB := ifelse(is.na(SYMBOL), GENEID, SYMBOL)]
setkey(features_db, 'GENEID')

# Save tables ----
# Goal: persist all objects needed by downstream scripts (DGEList, sample
# metadata, gene annotation table), saved both as RDS for R and as CSV for
# manual inspection.

# Save gene counts.
saveRDS(
  dge,
  file.path(output_path, paste0('DGEList_obj_', mapper,'.rds')))

write.csv(
  dge$counts,
  file = file.path(output_path, 'raw_matrix.csv'))

# Save metadata.
saveRDS(
  metadata,
  file.path(output_path, paste0('metadata_', mapper,'.rds')))

write.csv(
  metadata,
  file = file.path(output_path, paste0('metadata_', mapper,'.csv')),
  row.names = FALSE)

# Save gene information.
saveRDS(
  features_db,
  file.path(output_path, 'genesDb.rds'))

write.csv(
  features_db,
  file = file.path(output_path, 'genesDb.csv'),
  row.names = FALSE)

# Session info ----
sessioninfo::session_info()
