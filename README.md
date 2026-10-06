# tensorQTL Trans-QTL Workflow

This repository contains a [WDL](https://openwdl.org/) workflow for running trans-QTL mapping using [tensorQTL](https://github.com/broadinstitute/tensorqtl). It is designed to run on the [Terra](https://app.terra.bio/) platform and leverages GPU acceleration for large-scale trans-QTL analyses.

## Overview

Trans-QTL mapping tests associations between genetic variants and phenotypes (e.g., gene expression) across the entire genome, rather than only in the local genomic region of each phenotype. This workflow wraps the `tensorqtl` Python package in a cloud-ready WDL task, making it easy to run large-scale trans-QTL analyses on Terra or any Cromwell-compatible platform.

## Workflow

### `tensorqtl_trans_workflow` (`tensorQTL_trans.wdl`)

The workflow consists of a single task, `tensorqtl_trans`, which runs trans-QTL mapping using tensorQTL.

**Task: `tensorqtl_trans`**

Runs `python3 -m tensorqtl` in `--mode trans`, which tests all variant–phenotype pairs genome-wide and outputs nominal association statistics.

#### Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `plink_pgen` | File | PLINK2 `.pgen` genotype file |
| `plink_pvar` | File | PLINK2 `.pvar` variant information file |
| `plink_psam` | File | PLINK2 `.psam` sample information file |
| `phenotype_bed` | File | Phenotype file in compressed BED format (`.bed.gz`) (no index input required by this task) |
| `covariates` | File | Covariates file (tab-separated, samples as columns) |
| `interaction_file` | File? | Optional headerless TSV: sample ID and one numeric interaction value |
| `pval_threshold` | Float | Sparse-output threshold for nominal p-values (default `0.00001`); not an FDR cutoff |
| `batch_size` | Int | Variants per computation batch (default `1000`) |
| `prefix` | String | Output filename prefix |
| `maf_threshold` | Float | Minor allele frequency threshold for filtering variants |
| `fdr` | Float? | Legacy input; not applied in trans mode |
| `return_dense` | Boolean | If `true`, returns dense association matrices (not supported with interactions); if `false`, returns pairs below `pval_threshold` |
| `memory` | Int | Task memory in GB (default `120`) |
| `disk_space` | Int | Disk space to allocate (GB) |
| `num_threads` | Int | Number of CPU threads (default `32`) |
| `num_gpus` | Int | Number of NVIDIA L4 GPUs (default `1`; use `1` on this machine) |
| `num_preempt` | Int | Number of preemptible retries |

#### Outputs

| Output | Type | Description |
|--------|------|-------------|
| `trans_qtl` | File | Sparse nominal trans-QTL pairs in Parquet format (`<prefix>.trans_qtl_pairs.parquet`) |

#### Runtime

- **Docker image**: `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`
- **Machine**: `g2-standard-32`, with 32 vCPUs and 128 GB system RAM
- **GPU**: one NVIDIA L4 (`nvidia-l4`)
- **Task memory**: 120 GB by default, with RAM left for the operating system
- **GCP zone**: `us-central1-c`

When you update an existing Terra configuration, set `memory=120`,
`num_threads=32`, and `num_gpus=1`, or remove those values to use the new
defaults. Explicit input values override the defaults. The larger machine
increases system RAM. GPU memory remains 24 GB.

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

Phenotypes must be provided as a [gzipped](http://www.htslib.org/doc/bgzip.html) BED file (`.bed.gz`). The format expected by tensorQTL is:

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

1. Import `tensorQTL_trans.wdl` into Terra through a workflow repository or Dockstore.
2. Upload your input files to a Google Cloud Storage bucket.
3. Fill in the workflow inputs JSON with the GCS paths to your files and desired parameter values.
4. Launch the workflow on Terra.

## References

- [tensorQTL GitHub](https://github.com/broadinstitute/tensorqtl)
- [Taylor-Weiner et al., *Genome Biology* 2019](https://doi.org/10.1186/s13059-019-1851-8)
- [WDL specification](https://openwdl.org/)
- [PLINK2 documentation](https://www.cog-genomics.org/plink/2.0/)

## Optional interaction mapping

Set `tensorqtl_trans_workflow.tensorqtl_trans.interaction_file` to the GCS URI
of a two-column, headerless TSV. Leave it unset for ordinary trans mapping.
WDL preserves this input as `File?`, so Terra localizes it before validation.
The task checks that all required paths are readable. It also places the
three PLINK files under a common local prefix, even when Cromwell localizes
those files in separate directories.

Example CD4 interaction file (tab-separated):

```text
SAMPLE1	0.18
SAMPLE2	0.24
SAMPLE3	0.12
```

The file must contain one finite numeric value per listed sample. Sample IDs
must be nonempty and unique. Use a consistent fraction scale (for example,
0–1 for CD4 fractions).

When an interaction file is present, the task keeps only samples present in
all three files: the interaction TSV, the covariates TSV, and the phenotype
BED. It creates local copies with the same sample order as the original BED.
It removes samples outside this intersection from each copy. The four BED
metadata columns, phenotype rows, and covariate names remain in the copies.
The task streams the BED rows into a compressed BED file. It does not load
the complete phenotype matrix into memory for this step.

The task reports the number of retained samples and the number removed from
each file. It stops if there are no shared samples, or if the retained
interaction values do not vary. It checks for a duplicate interaction main
effect in the retained covariates. Retained samples must also be present in
the genotype files. Without an interaction file, the task uses the original
BED and covariates files.

For interaction runs, the task uses a small Python launcher for the tensorQTL
CLI. The launcher reads sample IDs in the aligned interaction TSV and the
genotype PSAM file as text.
This prevents numeric IDs from becoming numbers and preserves leading zeros.
It also preserves IDs such as `NA` as text. The reader change applies only to
those two files in that process. It does not change the installed
tensorQTL package or the sample IDs in any file.

The model includes genotype, the interaction variable's main effect, and
genotype × interaction variable, plus the covariates. For a CD4 scan:

- Put CD4 fractions in `interaction_file`.
- Put genetic PCs, expression PCs, technical covariates, and seven other
  fractions in `covariates`.
- Do not put CD4 fractions in `covariates`; the interaction option already
  includes their main effect. The task rejects a duplicate main effect,
  including a centered or rescaled copy.
- Omit one additional fraction because all nine fractions sum to one.
- Set `return_dense=false`. tensorQTL trans mode supports one interaction
  variable and sparse output only.

This tests whether a genetic effect varies with CD4 abundance. It does not
fit the joint nine-cell Decon-eQTL model or prove that an effect is exclusive
to CD4 cells.

The sparse output includes `pval_g`, `pval_i`, and `pval_gi`; use `pval_gi`
for the genotype × interaction test. The output threshold is a storage
filter, not multiple-testing correction. Apply a correction that accounts
for all tested trans pairs. With BED phenotype input, tensorQTL removes
pairs within ±5 Mb; dense ordinary output is not filtered this way.

`trans_qtl` is now optional: it is populated for sparse runs. For ordinary
dense runs, the workflow returns `trans_qtls_pval`, `trans_qtl_beta`,
`trans_qtl_beta_se`, and `trans_qtl_af` instead.

## Validation

GitHub Actions runs WDL validation, task command tests, a static check for
workflow-scope file-writing functions, and a real CPU tensorQTL 1.0.10 smoke
test. No Docker image is built. The command tests cover absent/present
interaction inputs, three-file sample intersection, plain and compressed
BED input, sample order, invalid values, safe path quoting,
unresolved cloud URIs, and separate PLINK localization directories.
The CPU smoke test compares an interaction run with partial sample overlap
against a manually filtered reference. It also checks ordinary sparse and
dense runs. The smoke inputs use numeric sample IDs with leading zeros.

The complete workflow has **not been tested on Terra**. The command tests
simulate cloud-to-local paths; they do not exercise Terra's localization
service. The CPU smoke uses a pinned tensorQTL package and does not confirm
the package version in the existing `:latest` GPU image.
