version 1.0

task tensorqtl_trans {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        File? interaction_file
        String prefix
        Float maf_threshold
        # Retained for compatibility; trans mode reports nominal p-values, not FDR.
        Float? fdr
        Boolean return_dense
        Float pval_threshold = 0.00001
        Int batch_size = 1000
        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_trans $*" >&2; }
        log "Validate localized inputs."
        plink_pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        plink_pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        plink_psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        interaction_file='~{sub(select_first([interaction_file, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$plink_pgen" "$plink_pvar" "$plink_psam" "$phenotype_bed" "$covariates"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        [[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { log "Input error: prefix must be a filename prefix."; exit 1; }
        if [[ -n "$interaction_file" ]]; then
            case "$interaction_file" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved interaction URI $interaction_file"; exit 1;; esac
            [[ -r "$interaction_file" ]] || { log "Input localization error: unreadable interaction file."; exit 1; }
            if ~{return_dense}; then
                log "Input error: tensorQTL trans interactions require sparse output (return_dense=false)."
                exit 1
            fi
            log "Clean samples: intersect interaction, covariates, and phenotype BED."
            python3 - "$interaction_file" "$phenotype_bed" "$covariates" <<'PY'
        import csv
        from datetime import datetime, timezone
        import gzip
        import math
        import sys

        def fail(message):
            sys.exit("Input error: " + message)

        interaction_path, phenotype_path, covariate_path = sys.argv[1:]
        opener = gzip.open if phenotype_path.endswith('.gz') else open
        with opener(phenotype_path, 'rt') as handle:
            header = next(csv.reader(handle, delimiter='\t'), [])
        bed_samples = header[4:]
        if not bed_samples or not all(bed_samples) or len(set(bed_samples)) != len(bed_samples):
            fail('phenotype BED must have nonempty, unique sample IDs after its four metadata columns.')
        values = {}
        with open(interaction_path) as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) != 2:
                    fail('interaction file must have exactly two columns: sample ID and value, without a header.')
                sample, text = row
                if not sample or sample in values:
                    fail('interaction sample IDs must be nonempty and unique.')
                try:
                    value = float(text)
                except ValueError:
                    fail('interaction values must be numeric; remove any header.')
                if not math.isfinite(value):
                    fail('interaction values must be finite.')
                values[sample] = value
        with open(covariate_path) as handle:
            cov_header = next(csv.reader(handle, delimiter='\t'), [])
        cov_samples = cov_header[1:]
        if not cov_samples or not all(cov_samples) or len(set(cov_samples)) != len(cov_samples):
            fail('covariate header must have nonempty, unique sample IDs after its ID column.')
        shared = set(values).intersection(cov_samples)
        samples = [sample for sample in bed_samples if sample in shared]
        if not samples:
            fail('no shared samples between interaction, covariates, and phenotype BED.')
        ordered = [values[sample] for sample in samples]
        if len(set(ordered)) < 2:
            fail('interaction values must vary across retained samples.')
        bed_index = {sample: index + 4 for index, sample in enumerate(bed_samples)}
        cov_index = {sample: index + 1 for index, sample in enumerate(cov_samples)}
        bed_columns = list(range(4)) + [bed_index[sample] for sample in samples]
        cov_columns = [0] + [cov_index[sample] for sample in samples]
        with opener(phenotype_path, 'rt') as handle, gzip.open('phenotype.aligned.bed.gz', 'wt') as output:
            reader = csv.reader(handle, delimiter='\t')
            writer = csv.writer(output, delimiter='\t', lineterminator='\n')
            writer.writerow([header[index] for index in bed_columns])
            next(reader)
            for row in reader:
                if len(row) != len(header):
                    fail('phenotype BED rows must match the header width.')
                writer.writerow([row[index] for index in bed_columns])
        interaction_mean = sum(ordered)/len(ordered)
        center = [value - interaction_mean for value in ordered]
        with open(covariate_path) as handle, open('covariates.aligned.tsv', 'w') as output:
            reader = csv.reader(handle, delimiter='\t')
            writer = csv.writer(output, delimiter='\t', lineterminator='\n')
            writer.writerow([cov_header[index] for index in cov_columns])
            next(reader)
            for row in reader:
                if len(row) != len(cov_header):
                    fail('covariate rows must match the header width.')
                aligned = [row[index] for index in cov_columns]
                try:
                    cov = [float(value) for value in aligned[1:]]
                except ValueError:
                    fail('covariates must be numeric.')
                if not all(math.isfinite(value) for value in cov):
                    fail('covariates must be finite.')
                # Detect the same main effect even if centered or expressed in percent.
                cov_mean = sum(cov)/len(cov)
                ccenter = [value - cov_mean for value in cov]
                scale = sum(a*b for a,b in zip(center, ccenter)) / sum(a*a for a in center)
                residual = sum((b-scale*a)**2 for a,b in zip(center, ccenter))
                if scale != 0 and residual <= 1e-12 * sum(b*b for b in ccenter):
                    fail('interaction main effect is already in the covariates: ' + row[0])
                writer.writerow(aligned)
        with open('interaction.aligned.tsv', 'w') as handle:
            for sample in samples:
                handle.write(sample + '\t' + str(values[sample]) + '\n')
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        print(f'[{timestamp}] stage=tensorqtl_trans Sample intersection: retained={len(samples)}; '
              f'BED dropped={len(bed_samples)-len(samples)}; '
              f'covariates dropped={len(cov_samples)-len(samples)}; '
              f'interaction dropped={len(values)-len(samples)}.', file=sys.stderr)
        PY
            phenotype_bed=phenotype.aligned.bed.gz
            covariates=covariates.aligned.tsv
        fi
        # Cromwell may localize PLINK components into separate directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode trans
              --covariates "$covariates" --maf_threshold ~{maf_threshold}
              --pval_threshold ~{pval_threshold} --batch_size ~{batch_size})
        if ~{return_dense}; then args+=(--return_dense); fi
        if [[ -n "$interaction_file" ]]; then args+=(--interaction interaction.aligned.tsv); fi
        if ~{defined(fdr)}; then log "Note: fdr is not applied in trans mode; saved p-values are nominal."; fi
        log "Start mapping; interaction enabled: ~{defined(interaction_file)}."
        if [[ -n "$interaction_file" ]]; then
            log "Read interaction sample IDs as text; preserve leading zeros."
            python3 - "${args[@]}" <<'PY'
        import os
        import runpy
        import sys
        import pandas as pd

        interaction_path = sys.argv[sys.argv.index('--interaction') + 1]
        original_read_csv = pd.read_csv

        def read_csv_with_sample_ids(path, *args, **kwargs):
            # tensorQTL's CLI infers interaction index types, unlike sample column names.
            # Change only the interaction file reader in this process; keep all IDs intact.
            if isinstance(path, (str, os.PathLike)) and os.fspath(path) == interaction_path:
                kwargs['dtype'] = {0: str}
                kwargs['keep_default_na'] = False
            return original_read_csv(path, *args, **kwargs)

        pd.read_csv = read_csv_with_sample_ids
        sys.argv[0] = 'tensorqtl'
        runpy.run_module('tensorqtl', run_name='__main__')
        PY
        else
            python3 -m tensorqtl "${args[@]}"
        fi
        log "Mapping completed."
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: num_threads
        preemptible: num_preempt
        predefinedMachineType: "g2-standard-16"
        gpuType: "nvidia-l4"
        gpuCount: num_gpus
        zones: ["us-central1-c"]
    }

    output {
        File? trans_qtl = prefix + ".trans_qtl_pairs.parquet"
        File? trans_qtls_pval = prefix + ".trans_qtl_pval.parquet"
        File? trans_qtl_beta = prefix + ".trans_qtl_beta.parquet"
        File? trans_qtl_beta_se = prefix + ".trans_qtl_beta_se.parquet"
        File? trans_qtl_af = prefix + ".trans_qtl_af.parquet"
    }
    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_trans_workflow {
    call tensorqtl_trans
    output {
        File? trans_qtl = tensorqtl_trans.trans_qtl
        File? trans_qtls_pval = tensorqtl_trans.trans_qtls_pval
        File? trans_qtl_beta = tensorqtl_trans.trans_qtl_beta
        File? trans_qtl_beta_se = tensorqtl_trans.trans_qtl_beta_se
        File? trans_qtl_af = tensorqtl_trans.trans_qtl_af
    }
}
