FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL org.opencontainers.image.description="R/Seurat/SCENT/bedtools container for single-cell enhancer-gene (ATAC x RNA) mapping with the SCENT package"
LABEL org.opencontainers.image.source="https://github.com/ccrobertson/scent-container"

# bedtools is a hard SystemRequirement of SCENT (used by CreatePeakToGeneList
# to intersect ATAC peaks against gene-body +/-500kb windows).
RUN apt-get update \
    && apt-get install -y --no-install-recommends bedtools procps \
    && rm -rf /var/lib/apt/lists/*

# CRAN deps: SCENT's own Imports, plus Seurat/optparse for the ETL wrapper
# scripts (SCENT itself never touches Seurat directly -- Seurat is only
# needed here to extract raw count matrices out of Seurat objects upstream
# of SCENT's own API).
RUN Rscript -e 'BiocManager::install(c( \
      "Seurat", "Matrix", "data.table", "dplyr", "stringr", "Hmisc", \
      "R.utils", "lme4", "boot", "optparse", "remotes", \
      "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", \
      "GenomicFeatures", "GenomicRanges", "IRanges", "GenomeInfoDb" \
    ), update = FALSE, ask = FALSE)'

RUN Rscript -e 'remotes::install_github("immunogenomics/SCENT@v1.0.1", upgrade = "never")'

COPY scripts/ /opt/scent/scripts/

# Bake a hg38 gene-body +/-500kb BED file into the image at build time, so
# CreatePeakToGeneList never needs network access at runtime on a compute
# node. See scripts/build_gene_bed.R for how it's derived.
RUN mkdir -p /opt/scent/ref \
    && Rscript /opt/scent/scripts/build_gene_bed.R /opt/scent/ref/hg38_genebody_500kb.bed \
    && test -s /opt/scent/ref/hg38_genebody_500kb.bed

WORKDIR /opt/scent
CMD ["R"]
