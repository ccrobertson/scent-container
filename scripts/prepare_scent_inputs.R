#!/usr/bin/env Rscript
#
# Builds a SCENT_obj.rds from this project's two separate Seurat objects
# (RNA and ATAC live in different objects, linked per-cell through the
# rna_barcode/atac_barcode columns in the metadata table -- 10x Multiome
# GEX and ATAC barcodes are NOT the same string, they're translated).
#
# Output: <outdir>/SCENT_obj.rds, ready to feed to run_scent_chunk.R, plus
# a small chunk_sizes.txt summary of how many gene-peak pairs landed in
# each parallelization chunk.

suppressMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(SCENT)
})

opt_list <- list(
  make_option("--rna_rds", type = "character",
              help = "Seurat object (RNA), raw counts assay/layer"),
  make_option("--atac_rds", type = "character",
              help = "Seurat object (ATAC peaks x cells), raw counts assay/layer"),
  make_option("--meta_csv", type = "character", default = "",
              help = "Metadata CSV with a row per cell [default: use rna_rds's own meta.data]"),
  make_option("--rna_assay", type = "character", default = "RNA"),
  make_option("--rna_layer", type = "character", default = "counts"),
  make_option("--atac_assay", type = "character", default = "RNA"),
  make_option("--atac_layer", type = "character", default = "counts"),
  make_option("--cell_col_rna", type = "character", default = "rna_barcode",
              help = "Column in metadata giving the RNA object's colname for each cell"),
  make_option("--cell_col_atac", type = "character", default = "atac_barcode",
              help = "Column in metadata giving the ATAC object's colname for each cell"),
  make_option("--qc_col", type = "character", default = "pass_all_filters",
              help = "Metadata column to require TRUE; set to '' to skip"),
  make_option("--singlet_col", type = "character", default = "droplet_type_strict",
              help = "Metadata column to filter on; set to '' to skip"),
  make_option("--singlet_value", type = "character", default = "SNG"),
  make_option("--genebed", type = "character",
              default = "/opt/scent/ref/hg38_genebody_500kb.bed",
              help = "Gene-body +/-500kb BED file for CreatePeakToGeneList"),
  make_option("--nbatch", type = "integer", default = 200,
              help = "Number of gene-peak chunks to split into for array parallelization"),
  make_option("--outdir", type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$rna_rds), !is.null(opt$atac_rds), !is.null(opt$outdir))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  as.logical(toupper(as.character(x)))
}

cat("Loading RNA object:", opt$rna_rds, "\n")
rna_obj <- readRDS(opt$rna_rds)
rna_counts <- GetAssayData(rna_obj, assay = opt$rna_assay, layer = opt$rna_layer)

cat("Loading ATAC object:", opt$atac_rds, "\n")
atac_obj <- readRDS(opt$atac_rds)
atac_counts <- GetAssayData(atac_obj, assay = opt$atac_assay, layer = opt$atac_layer)

if (nzchar(opt$meta_csv)) {
  cat("Loading metadata:", opt$meta_csv, "\n")
  meta <- as.data.frame(fread(opt$meta_csv))
} else {
  cat("Using rna_rds's own meta.data (no --meta_csv given)\n")
  meta <- rna_obj@meta.data
  meta[[opt$cell_col_rna]] <- rownames(meta)
}

for (col in c(opt$cell_col_rna, opt$cell_col_atac)) {
  if (!col %in% colnames(meta)) stop("Column '", col, "' not found in metadata")
}

n0 <- nrow(meta)
if (nzchar(opt$qc_col)) {
  stopifnot(opt$qc_col %in% colnames(meta))
  meta <- meta[to_logical(meta[[opt$qc_col]]) %in% TRUE, ]
  cat("After", opt$qc_col, "filter:", nrow(meta), "/", n0, "cells\n")
}
if (nzchar(opt$singlet_col)) {
  stopifnot(opt$singlet_col %in% colnames(meta))
  meta <- meta[as.character(meta[[opt$singlet_col]]) == opt$singlet_value, ]
  cat("After", opt$singlet_col, "==", opt$singlet_value, "filter:", nrow(meta), "cells\n")
}

keep <- meta[[opt$cell_col_rna]] %in% colnames(rna_counts) &
        meta[[opt$cell_col_atac]] %in% colnames(atac_counts)
cat("Cells with matches in both RNA and ATAC matrices:", sum(keep), "/", nrow(meta), "\n")
meta <- meta[keep, ]
meta <- meta[!duplicated(meta[[opt$cell_col_rna]]), ]

rna_counts <- rna_counts[, meta[[opt$cell_col_rna]]]
atac_counts <- atac_counts[, meta[[opt$cell_col_atac]]]
colnames(atac_counts) <- meta[[opt$cell_col_rna]]  # canonical cell id = rna_barcode

rna_counts <- as(as(rna_counts, "CsparseMatrix"), "dgCMatrix")
atac_counts <- as(as(atac_counts, "CsparseMatrix"), "dgCMatrix")

meta$cell <- meta[[opt$cell_col_rna]]
rownames(meta) <- meta$cell
if ("rna_umis" %in% colnames(meta)) {
  meta$log_nUMI <- log(as.numeric(meta$rna_umis))
} else {
  warning("No 'rna_umis' column found -- log_nUMI covariate will not be available")
}

cat("Final matched cell count:", nrow(meta), "\n")
cat("RNA counts:", nrow(rna_counts), "genes x", ncol(rna_counts), "cells\n")
cat("ATAC counts:", nrow(atac_counts), "peaks x", ncol(atac_counts), "cells\n")

cat("Constructing SCENT object...\n")
scent_obj <- CreateSCENTObj(
  rna = rna_counts, atac = atac_counts, meta.data = meta,
  peak.info = data.frame(gene = character(0), peak = character(0)),
  covariates = character(0),  # covariates chosen at run_scent_chunk.R time
  celltypes = "cell_type"
)

cat("Building gene-peak pair chunks via CreatePeakToGeneList (bedtools intersect)...\n")
scent_obj <- CreatePeakToGeneList(
  scent_obj, genebed = opt$genebed, nbatch = opt$nbatch,
  tmpfile = file.path(opt$outdir, "tmp_atac_peak.bed"),
  intersectedfile = file.path(opt$outdir, "tmp_atac_peak_intersected.bed.gz")
)

out_rds <- file.path(opt$outdir, "SCENT_obj.rds")
saveRDS(scent_obj, out_rds)

chunk_sizes <- vapply(scent_obj@peak.info.list, nrow, integer(1))
writeLines(
  paste0("chunk_", seq_along(chunk_sizes), "\t", chunk_sizes, "\tpairs"),
  file.path(opt$outdir, "chunk_sizes.txt")
)

cat("Wrote", out_rds, "with", length(chunk_sizes), "chunks",
    "(", sum(chunk_sizes), "total gene-peak pairs )\n")
cat("Cell type counts:\n")
print(table(meta$cell_type))
