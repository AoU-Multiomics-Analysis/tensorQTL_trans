# Chromosome trans analysis implementation plan

**Goal:** Split the input PLINK files by chromosome, run trans mapping for each chromosome, and merge the results.

**Architecture:** A CPU preparation task performs the existing sample intersection once. A CPU PLINK2 task creates chromosome files. A WDL scatter maps each chromosome against all retained phenotypes. A CPU task streams each result into the final Parquet files.

**Constraints:** WDL 1.0; typed File inputs until localization; task-only file creation; safe CLI arguments; command logs; no cloud submission; no local Docker build.

**Review focus:** Empty chromosome results; numeric sample IDs; exact chromosome labels and variant order; dense matrix indexes; localized files with spaces and cloud URI rejection.

## Tasks

- [ ] Add command tests for preparation, chromosome split, scatter dependencies, merge, and logs. Verify the new tests fail before implementation.
- [ ] Move sample cleaning into one task. Split all chromosomes present in PVAR with pinned PLINK2, preserving samples and variant metadata. Keep every file typed in the scatter.
- [ ] Stream sparse and dense results into the existing final output names. Preserve empty sparse results and fail on inconsistent schemas.
- [ ] Extend the GitHub Actions CPU smoke test to compare merged chromosome results with whole-genome results, including empty results and dense output. Run WDL and command checks locally.
- [ ] Obtain an independent code review, address material findings, and open a pull request. Wait for GitHub Actions results.

**Validation limit:** The full workflow must be tested on Terra separately. This change reduces genotype memory per GPU task; each task still loads all phenotypes.
