# tensorQTL Trans-QTL Workflow

This repository contains a [WDL](https://openwdl.org/) workflow for running trans-QTL mapping using [tensorQTL](https://github.com/broadinstitute/tensorqtl). It is designed to run on the [Terra](https://app.terra.bio/) platform and leverages GPU acceleration for large-scale trans-QTL analyses.

## Overview

Trans-QTL mapping tests associations between genetic variants and phenotypes (e.g., gene expression) across the entire genome, rather than only in the local genomic region of each phenotype. This workflow wraps the `tensorqtl` Python package in a cloud-ready WDL task, making it easy to run large-scale trans-QTL analyses on Terra or any Cromwell-compatible platform.

## Workflow

### `tensorqtl_trans_workflow` (`tensorQTL_trans.wdl`)

The workflow consists of a single task, `tensorqtl_trans`, which runs trans-QTL mapping using tensorQTL.

**Task: `tensorqtl_trans`**

Runs `python3 -m tensorqtl` in `--mode trans`, which tests all variant–phenotype pairs genome-wide and outputs significant associations.

#### Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `plink_pgen` | File | PLINK2 `.pgen` genotype file |
| `plink_pvar` | File | PLINK2 `.pvar` variant information file |
| `plink_psam` | File | PLINK2 `.psam` sample information file |
| `phenotype_bed` | File | Phenotype file in compressed BED format (`.bed.gz`) with corresponding index (`.bed.gz.tbi`) |
| `covariates` | File | Covariates file (tab-separated, samples as columns) |
| `prefix` | String | Output filename prefix |
| `maf_threshold` | Float | Minor allele frequency threshold for filtering variants |
| `fdr` | Float? | (Optional) FDR threshold for reporting significant associations |
| `return_dense` | Boolean | If `true`, returns all variant–phenotype pairs; if `false`, returns only significant pairs |
| `memory` | Int | Memory to allocate (GB) |
| `disk_space` | Int | Disk space to allocate (GB) |
| `num_threads` | Int | Number of CPU threads |
| `num_gpus` | Int | Number of GPUs (NVIDIA Tesla P100) |
| `num_preempt` | Int | Number of preemptible retries |

#### Outputs

| Output | Type | Description |
|--------|------|-------------|
| `trans_qtl` | File | Significant trans-QTL pairs in Parquet format (`<prefix>.trans_qtl_pairs.parquet`) |

#### Runtime

- **Docker image**: `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`
- **GPU**: NVIDIA Tesla P100 (`nvidia-tesla-p100`)
- **GCP zone**: `us-central1-c`

## Data Preparation

### Genotype Data (PLINK2 format)

Genotype data must be provided in [PLINK2](https://www.cog-genomics.org/plink/2.0/) binary format (`.pgen`/`.pvar`/`.psam`). To convert from VCF:

```bash
plink2 --vcf input.vcf.gz \
       --make-pgen \
       --out output_prefix \
       --max-alleles 2 \
       --maf 0.01
```

Recommended preprocessing steps:
- Retain only biallelic SNPs (`--max-alleles 2`)
- Apply an initial MAF filter (the workflow also applies a MAF filter at runtime via `maf_threshold`)
- Ensure sample IDs in `.psam` match those in the phenotype and covariates files

### Phenotype Data (BED format)

Phenotypes must be provided as a [bgzipped](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html)-indexed BED file (`.bed.gz` + `.bed.gz.tbi`). The format expected by tensorQTL is:

- **Tab-separated**
- First four columns: `#chr`, `start`, `end`, `phenotype_id`
- Remaining columns: one per sample (sample IDs must match genotype data)
- Rows represent individual phenotypes (e.g., genes)

To prepare:

```bash
# Sort and bgzip
sort -k1,1 -k2,2n phenotypes.bed | bgzip -c > phenotypes.bed.gz

# Index with tabix
tabix -p bed phenotypes.bed.gz
```

### Covariates

The covariates file should be tab-separated with:
- First column: covariate name
- Remaining columns: one per sample (sample IDs must match genotype data)

Common covariates include:
- Genotype principal components (PCs)
- Phenotype PCs (e.g., PEER factors or expression PCs)
- Technical covariates (e.g., sequencing batch, sex, age)

Example format:

```
ID    SAMPLE1    SAMPLE2    SAMPLE3
PC1   0.012      -0.034     0.021
PC2   0.005      0.011      -0.009
sex   1          2          1
```

### Sample Consistency

Before running the workflow, ensure that sample IDs are consistent across all three input files:
- PLINK2 `.psam` (second column: `IID`)
- Phenotype BED header (columns 5 onward)
- Covariates header (columns 2 onward)

## Running on Terra

1. Import the workflow from [Dockstore](https://dockstore.org/) using the `.dockstore.yml` configuration.
2. Upload your input files to a Google Cloud Storage bucket.
3. Fill in the workflow inputs JSON with the GCS paths to your files and desired parameter values.
4. Launch the workflow on Terra.

## References

- [tensorQTL GitHub](https://github.com/broadinstitute/tensorqtl)
- [Taylor-Weiner et al., *Genome Biology* 2019](https://doi.org/10.1186/s13059-019-1851-8)
- [WDL specification](https://openwdl.org/)
- [PLINK2 documentation](https://www.cog-genomics.org/plink/2.0/)
