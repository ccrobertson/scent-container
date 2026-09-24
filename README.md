# scent-container

A container image for running [SCENT](https://github.com/immunogenomics/SCENT)
(Single-Cell ENhancer Target gene mapping) on single-cell multiome (RNA + ATAC)
data, including from Seurat objects.

Built on `bioconductor/bioconductor_docker` with:

- R packages: `SCENT` (v1.0.1, from GitHub), `Seurat`, and SCENT's own
  dependencies (`data.table`, `lme4`, `stringr`, `boot`, `MASS`, `Matrix`,
  `Hmisc`, `R.utils`)
- `bedtools` (a hard `SystemRequirement` of SCENT, used to build cis
  gene-peak windows)
- A pre-built hg38 gene-body ±500kb BED file baked in at
  `/opt/scent/ref/hg38_genebody_500kb.bed` (see `scripts/build_gene_bed.R`),
  so no network access is needed at analysis runtime.

## Build

Images are built and pushed automatically to
[ghcr.io/ccrobertson/scent-container](https://github.com/ccrobertson/scent-container/pkgs/container/scent-container)
by GitHub Actions (`.github/workflows/build.yml`) on every push to `main`
that touches `Dockerfile` or `scripts/`.

## Use on an HPC cluster with Singularity/Apptainer

No local build/root needed -- pull the image GitHub Actions already built:

```bash
module load singularity   # or apptainer, depending on your cluster
singularity pull scent.sif docker://ghcr.io/ccrobertson/scent-container:latest
```

**Always run with `--no-home`.** Singularity mounts your host `$HOME` by
default, and if you have a personal R library there (e.g. from prior
`module load R/...` work, or an `.Rprofile` that customizes `.libPaths()`),
R will find and load *those* host packages ahead of the container's own --
concretely, this breaks `library(Seurat)` if your host has an older
`ggplot2` than Seurat's container version requires. `--no-home` prevents
the host library from ever shadowing the container's:

```bash
singularity exec --no-home scent.sif Rscript -e 'library(SCENT); library(Seurat)'
```

Paths outside `$HOME` (e.g. `/nfs/turbo/...`) are unaffected and still
reachable with `--no-home`, as long as your cluster's Singularity config
auto-binds them (Great Lakes does for `/nfs`).

**BLAS/OpenMP threads are pinned to 1** (`OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `BLAS_NUM_THREADS`, set as image
`ENV` vars). This matters if you use SCENT's `ncores` argument (or any
`parallel::mclapply`/`boot(parallel="multicore")` code) for outer-level
parallelism: without pinning, each forked process can *also* multithread
its own linear algebra internally, oversubscribing the CPUs your job
actually has allocated (e.g. 4 forked processes x 2 internal threads each
= 8 threads contending for 4 allocated cores) -- this silently cancels out
the parallel speedup with no error or warning; `ncores=4` ends up no
faster than `ncores=1`. Confirmed empirically on a SCENT bootstrap
workload: unpinned, 4-way parallel took the same wall time as serial;
pinned, it was ~2.8x faster. If you override these env vars for some other
reason, be aware you may reintroduce this.

Note: `OMP_NUM_THREADS`/`OMP_THREAD_LIMIT` specifically are also patched
directly in `/usr/local/lib/R/etc/Renviron.site` (not just set via Docker
`ENV`), because the base `bioconductor_docker` image's own `Renviron.site`
(for Bioconductor's build-machine testing) hardcodes
`OMP_NUM_THREADS=2`/`OMP_THREAD_LIMIT=2`, and R's `Renviron.site`
processing overwrites the process environment at startup -- a Docker `ENV`
alone is silently undone from R's perspective. `OPENBLAS_NUM_THREADS`/
`MKL_NUM_THREADS`/`BLAS_NUM_THREADS` aren't touched by that file, so `ENV`
was already sufficient for those.

## Scripts

- `scripts/build_gene_bed.R` -- generates the hg38 gene-body ±500kb BED file
  baked into the image at build time.
- `scripts/prepare_scent_inputs.R` -- ETL helper: given an RNA Seurat object,
  an ATAC (peaks x cells) Seurat object, and a per-cell metadata table, builds
  a `SCENT_obj.rds` ready for `SCENT_algorithm()`, via SCENT's own
  `CreateSCENTObj()`/`CreatePeakToGeneList()`. Run it with:

  ```bash
  singularity exec --no-home scent.sif Rscript /opt/scent/scripts/prepare_scent_inputs.R \
    --rna_rds path/to/rna.rds \
    --atac_rds path/to/atac.rds \
    --meta_csv path/to/metadata.csv \
    --outdir path/to/output
  ```

  Run `Rscript /opt/scent/scripts/prepare_scent_inputs.R --help` inside the
  container for the full list of options (assay/layer names, barcode-matching
  columns, QC filters, genome BED path, number of parallelization chunks).
