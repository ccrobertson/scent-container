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
[ghcr.io/ccrobertson/scent](https://github.com/ccrobertson/scent-container/pkgs/container/scent)
by GitHub Actions (`.github/workflows/build.yml`) on every push to `main`
that touches `Dockerfile` or `scripts/`.

## Use on an HPC cluster with Singularity/Apptainer

No local build/root needed -- pull the image GitHub Actions already built:

```bash
module load singularity   # or apptainer, depending on your cluster
singularity pull scent.sif docker://ghcr.io/ccrobertson/scent:latest
```

## Scripts

- `scripts/build_gene_bed.R` -- generates the hg38 gene-body ±500kb BED file
  baked into the image at build time.
- `scripts/prepare_scent_inputs.R` -- ETL helper: given an RNA Seurat object,
  an ATAC (peaks x cells) Seurat object, and a per-cell metadata table, builds
  a `SCENT_obj.rds` ready for `SCENT_algorithm()`, via SCENT's own
  `CreateSCENTObj()`/`CreatePeakToGeneList()`. Run it with:

  ```bash
  singularity exec scent.sif Rscript /opt/scent/scripts/prepare_scent_inputs.R \
    --rna_rds path/to/rna.rds \
    --atac_rds path/to/atac.rds \
    --meta_csv path/to/metadata.csv \
    --outdir path/to/output
  ```

  Run `Rscript /opt/scent/scripts/prepare_scent_inputs.R --help` inside the
  container for the full list of options (assay/layer names, barcode-matching
  columns, QC filters, genome BED path, number of parallelization chunks).
