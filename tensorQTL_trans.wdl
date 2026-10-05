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
            python3 - "$interaction_file" "$phenotype_bed" "$covariates" <<'PY'
        import csv
        import gzip
        import math
        import sys

        def fail(message):
            sys.exit("Input error: " + message)

        interaction_path, phenotype_path, covariate_path = sys.argv[1:]
        opener = gzip.open if phenotype_path.endswith('.gz') else open
        with opener(phenotype_path, 'rt') as handle:
            header = next(csv.reader(handle, delimiter='\t'))
        samples = header[4:]
        if not samples or len(set(samples)) != len(samples):
            fail('phenotype BED must have unique sample IDs after its four metadata columns.')
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
        if set(values) != set(samples):
            fail('interaction sample IDs must match the phenotype BED exactly.')
        ordered = [values[sample] for sample in samples]
        if len(set(ordered)) < 2:
            fail('interaction values must vary across samples.')
        with open(covariate_path) as handle:
            reader = csv.reader(handle, delimiter='\t')
            cov_header = next(reader)
            if cov_header[1:] != samples:
                fail('covariate sample IDs and order must match the phenotype BED.')
            for row in reader:
                if len(row) != len(samples) + 1:
                    fail('covariate rows must match the header width.')
                try:
                    cov = [float(value) for value in row[1:]]
                except ValueError:
                    fail('covariates must be numeric.')
                if not all(math.isfinite(value) for value in cov):
                    fail('covariates must be finite.')
                # Detect the same main effect even if centered or expressed in percent.
                center = [value - sum(ordered)/len(ordered) for value in ordered]
                ccenter = [value - sum(cov)/len(cov) for value in cov]
                scale = sum(a*b for a,b in zip(center, ccenter)) / sum(a*a for a in center)
                residual = sum((b-scale*a)**2 for a,b in zip(center, ccenter))
                if scale != 0 and residual <= 1e-12 * sum(b*b for b in ccenter):
                    fail('interaction main effect is already in the covariates: ' + row[0])
        with open('interaction.aligned.tsv', 'w') as handle:
            for sample in samples:
                handle.write(sample + '\t' + str(values[sample]) + '\n')
        PY
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
        python3 -m tensorqtl "${args[@]}"
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
