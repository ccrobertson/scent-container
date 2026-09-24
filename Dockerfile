FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL org.opencontainers.image.description="R/Seurat/SCENT/bedtools container for single-cell enhancer-gene (ATAC x RNA) mapping with the SCENT package"
LABEL org.opencontainers.image.source="https://github.com/ccrobertson/scent-container"

# SCENT's own bootstrap parallelizes via parallel::mclapply (forking
# `ncores` processes). Without pinning these, R's BLAS/LAPACK backend can
# ALSO multithread each glm.fit() internally (e.g. 2 OpenMP threads by
# default in this image), so N forked processes x M internal threads each
# oversubscribes the CPUs actually allocated (e.g. 4 forks x 2 threads = 8
# threads competing for 4 allocated cores) -- this silently cancels out
# the benefit of the outer-level fork parallelism, making `ncores=4` run no
# faster than `ncores=1` with no error or warning. Confirmed empirically:
# unpinned, 4-way parallel bootstrap took the same wall time as serial;
# pinned to 1 thread each, 4-way parallel was ~2.8x faster than serial.
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV BLAS_NUM_THREADS=1

# The base image's own /usr/local/lib/R/etc/Renviron.site (set up for
# Bioconductor's build-machine testing, unrelated to us) hardcodes
# OMP_NUM_THREADS=2 and OMP_THREAD_LIMIT=2. R's Renviron.site processing
# overwrites the process environment at startup, silently undoing the ENV
# line above from R's own perspective (confirmed: Sys.getenv("OMP_NUM_THREADS")
# read "2" even with ENV OMP_NUM_THREADS=1 set) -- so a Docker ENV alone
# isn't sufficient here. Patch it directly rather than relying on ENV.
RUN sed -i \
      -e 's/^OMP_NUM_THREADS=.*/OMP_NUM_THREADS=1/' \
      -e 's/^OMP_THREAD_LIMIT=.*/OMP_THREAD_LIMIT=1/' \
      /usr/local/lib/R/etc/Renviron.site

# bedtools is a hard SystemRequirement of SCENT (used by CreatePeakToGeneList
# to intersect ATAC peaks against gene-body +/-500kb windows).
RUN apt-get update \
    && apt-get install -y --no-install-recommends bedtools procps \
    && rm -rf /var/lib/apt/lists/*

# The base image ships CRAN packages pinned to whatever was current at the
# Bioconductor release date (e.g. ggplot2 3.5.1), which can be too old for
# versioned NAMESPACE imports declared by packages installed later (Seurat
# requires ggplot2 >= 3.5.2). install.packages() doesn't upgrade an
# already-installed dependency just because a new package's NAMESPACE
# requires a newer version -- it installs fine but fails at library() time.
# Refresh all preinstalled CRAN packages up front to avoid that whack-a-mole.
RUN Rscript -e 'update.packages(ask = FALSE, checkBuilt = TRUE, repos = BiocManager::repositories())'

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
