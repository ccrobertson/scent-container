#!/usr/bin/env Rscript
# Builds a BED4 file of hg38 gene bodies padded +/-500kb (clipped to
# chromosome bounds), for use as the `genebed` argument to SCENT's
# CreatePeakToGeneList(). Run once at Docker build time -- no internet
# access needed at analysis runtime.
#
# Usage: Rscript build_gene_bed.R <output.bed>

args <- commandArgs(trailingOnly = TRUE)
out_bed <- if (length(args) >= 1) args[1] else "hg38_genebody_500kb.bed"

suppressMessages({
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomeInfoDb)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
std_chroms <- paste0("chr", c(1:22, "X", "Y"))

g <- genes(txdb)
g <- g[as.character(GenomeInfoDb::seqnames(g)) %in% std_chroms]

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(g),
                              keytype = "ENTREZID", column = "SYMBOL",
                              multiVals = "first")
g$gene <- sym
g <- g[!is.na(g$gene)]

body_bed <- data.frame(
  chr   = as.character(GenomeInfoDb::seqnames(g)),
  start = pmax(start(g) - 1L, 0L),   # BED is 0-based half-open
  end   = end(g),
  gene  = g$gene
)

genome_file <- tempfile(fileext = ".genome")
body_file   <- tempfile(fileext = ".bed")
sorted_file <- tempfile(fileext = ".bed")

sl <- GenomeInfoDb::seqlengths(txdb)[std_chroms]
write.table(data.frame(names(sl), as.integer(sl)), genome_file,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(body_bed, body_file,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

system(paste("sort -k1,1 -k2,2n", body_file, ">", sorted_file))

# Pad +/-500kb, clipped to chromosome bounds via the .genome file.
status <- system(paste(
  "bedtools slop -i", sorted_file, "-g", genome_file, "-b 500000 >", out_bed
))
if (status != 0 || !file.exists(out_bed) || file.info(out_bed)$size == 0) {
  stop("Failed to build gene-body BED file at ", out_bed)
}

cat("Wrote", length(readLines(out_bed)), "gene-body +/-500kb windows to", out_bed, "\n")
