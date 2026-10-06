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
        Int memory = 120
        Int disk_space
        Int num_threads = 32
        Int num_gpus = 1
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
        # Cromwell may localize PLINK components into separate directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode trans
              --covariates "$covariates" --maf_threshold ~{maf_threshold}
              --pval_threshold ~{pval_threshold} --batch_size ~{batch_size})
        if ~{return_dense}; then args+=(--return_dense); fi
        if [[ -n "$interaction_file" ]]; then args+=(--interaction "$interaction_file"); fi
        if ~{defined(fdr)}; then log "Note: fdr is not applied in trans mode; saved p-values are nominal."; fi
        log "Resources: machine=g2-standard-32; memory=~{memory}GB; CPUs=~{num_threads}; GPUs=~{num_gpus}."
        if [[ -r /proc/meminfo ]]; then awk '/^MemTotal:/ {print "Host " $0 > "/dev/stderr"}' /proc/meminfo; fi
        for limit in /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory/memory.limit_in_bytes; do
            if [[ -r "$limit" ]]; then log "Container memory limit (bytes): $(cat "$limit")"; break; fi
        done
        log "Start mapping; interaction enabled: ~{defined(interaction_file)}."
        if [[ -n "$interaction_file" ]]; then
            log "Read interaction and genotype sample IDs as text; preserve leading zeros."
            python3 -u - "${args[@]}" 1>&2 <<'PY'
        import os
        import runpy
        import sys
        import pandas as pd

        interaction_path = sys.argv[sys.argv.index('--interaction') + 1]
        psam_path = sys.argv[1] + '.psam'
        original_read_csv = pd.read_csv

        def read_csv_with_sample_ids(path, *args, **kwargs):
            # tensorQTL's CLI infers interaction index types, unlike sample column names.
            # Change only the interaction and PSAM readers in this process; keep IDs intact.
            if isinstance(path, (str, os.PathLike)) and os.fspath(path) == interaction_path:
                kwargs['dtype'] = {0: str}
                kwargs['keep_default_na'] = False
            elif isinstance(path, (str, os.PathLike)) and os.fspath(path) == psam_path:
                dtype = kwargs.get('dtype')
                if dtype is None or isinstance(dtype, dict):
                    kwargs['dtype'] = {**(dtype or {}), 0: str, '#IID': str, 'IID': str,
                                       '#FID': str, 'FID': str}
                kwargs['keep_default_na'] = False
            return original_read_csv(path, *args, **kwargs)

        pd.read_csv = read_csv_with_sample_ids
        from tensorqtl import trans
        original_map_trans = trans.map_trans

        def map_trans_with_empty_results(*args, **kwargs):
            try:
                return original_map_trans(*args, **kwargs)
            except ValueError as error:
                # tensorQTL 1.0.10 concatenates empty hit lists when no interaction
                # pair passes the threshold. Handle only that specific empty result.
                frame = error.__traceback__
                while frame is not None:
                    state = frame.tb_frame.f_locals
                    if (frame.tb_frame.f_code is original_map_trans.__code__
                            and str(error) == 'need at least one array to concatenate'
                            and state.get('interaction_s') is not None
                            and all(name in state and len(state[name]) == 0 for name in
                                    ['tstat_g_list', 'tstat_i_list', 'tstat_gi_list', 'af_list', 'ix0', 'ix1'])):
                        print('No interaction pairs passed the output threshold.', file=sys.stderr)
                        return pd.DataFrame({
                            'variant_id': pd.Series(dtype='string'), 'phenotype_id': pd.Series(dtype='string'),
                            'pval_g': pd.Series(dtype='float64'), 'pval_i': pd.Series(dtype='float64'),
                            'pval_gi': pd.Series(dtype='float64'), 'af': pd.Series(dtype='float32')})
                    frame = frame.tb_next
                raise

        trans.map_trans = map_trans_with_empty_results
        sys.argv[0] = 'tensorqtl'
        runpy.run_module('tensorqtl', run_name='__main__')
        PY
        else
            python3 -u -m tensorqtl "${args[@]}" 1>&2
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
        predefinedMachineType: "g2-standard-32"
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

task prepare_samples {
    input {
        File phenotype_bed
        File covariates
        File? interaction_file
        String prefix
        Boolean return_dense
        Int disk_space
        Int auxiliary_memory = 16
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_trans $*" >&2; }
        log "Validate localized inputs."
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        interaction_file='~{sub(select_first([interaction_file, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$phenotype_bed" "$covariates"; do
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
        if [[ -n "$interaction_file" ]]; then cp interaction.aligned.tsv prepared.interaction.tsv; fi
        log "Sample preparation completed."
    >>>
    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "~{auxiliary_memory}GB"
        cpu: 2
        disks: "local-disk ~{disk_space} HDD"
        preemptible: num_preempt
    }
    output {
        # Return original files in ordinary mode so BED compression is detected correctly.
        File bed = if defined(interaction_file) then "phenotype.aligned.bed.gz" else phenotype_bed
        File cov = if defined(interaction_file) then "covariates.aligned.tsv" else covariates
        File? interaction = if defined(interaction_file) then "prepared.interaction.tsv" else interaction_file
    }
}

task split_chromosomes {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        Int disk_space
        Int split_memory = 8
        Int split_threads = 4
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=split_chromosomes $*" >&2; }
        pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        for path in "$pgen" "$pvar" "$psam"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        [[ ~{split_memory} -ge 2 && ~{split_threads} -ge 1 ]] || { log "Input error: split_memory must be >=2 and split_threads must be positive."; exit 1; }
        log "Discover chromosomes in PVAR."
        awk '!/^#/ && NF {if (!seen[$1]++) print $1}' "$pvar" > chromosomes.txt
        [[ -s chromosomes.txt ]] || { log "Input error: PVAR contains no variants."; exit 1; }
        # Compare sample identities as text, including their order and leading zeros.
        sample_ids() {
            awk 'NR==1 {for(i=1;i<=NF;i++) {if($i=="#IID" || $i=="IID") iid=i; if($i=="#FID" || $i=="FID") fid=i; if($i=="SID") sid=i} if(!iid) exit 1; next}
                 {print (fid ? $fid : "0") "\t" $iid "\t" (sid ? $sid : "0")}' "$1"
        }
        sample_ids "$psam" > input.samples
        mkdir chromosomes
        index=0
        while IFS= read -r chromosome; do
            [[ "$chromosome" =~ ^[A-Za-z0-9_.-]+$ && "$chromosome" != -* ]] || { log "Input error: unsupported chromosome label $chromosome"; exit 1; }
            index=$((index+1))
            output=$(printf 'chromosomes/chr_%06d' "$index")
            log "Split chromosome=$chromosome; output=$output."
            plink2 --pgen "$pgen" --pvar "$pvar" --psam "$psam" \
                --chr "$chromosome" --allow-extra-chr --make-pgen \
                --memory $((~{split_memory}*1000-1000)) --threads ~{split_threads} --out "$output" 1>&2
            sample_ids "$output.psam" > split.samples
            diff input.samples split.samples || { log "Split error: sample IDs or order changed."; exit 1; }
            # PLINK can normalize chromosome names. Restore the exact original metadata
            # only after verifying variant IDs, positions, alleles, and row order.
            awk -v chr="$chromosome" '/^#/ || $1==chr' "$pvar" > original.pvar
            awk '!/^#/ && NF {print $2,$3,$4,$5}' original.pvar > original.variants
            awk '!/^#/ && NF {print $2,$3,$4,$5}' "$output.pvar" > split.variants
            diff original.variants split.variants || { log "Split error: variants or order changed."; exit 1; }
            mv original.pvar "$output.pvar"
            cp "$psam" "$output.psam"
        done < chromosomes.txt
        log "Split completed: $index chromosomes."
    >>>
    runtime {
        docker: "quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
        memory: "~{split_memory}GB"
        cpu: split_threads
        disks: "local-disk ~{disk_space} HDD"
        preemptible: num_preempt
    }
    output {
        Array[String] chromosomes = read_lines("chromosomes.txt")
        Array[File] pgens = glob("chromosomes/chr_*.pgen")
        Array[File] pvars = glob("chromosomes/chr_*.pvar")
        Array[File] psams = glob("chromosomes/chr_*.psam")
    }
}

task merge_trans {
    input {
        Array[File] pairs
        Array[File] pvals
        Array[File] betas
        Array[File] beta_ses
        Array[File] afs
        Int chromosome_count
        String prefix
        Boolean return_dense
        Int disk_space
        Int auxiliary_memory = 16
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        log() { echo "[$(date -u +%FT%TZ)] stage=merge_trans $*" >&2; }
        log "Merge chromosome results."
        # Convert typed File arrays to local paths only during command rendering.
        printf '%s' '~{sub(sep("\n", pairs), "'", "'\"'\"'")}' > pairs.list
        printf '%s' '~{sub(sep("\n", pvals), "'", "'\"'\"'")}' > pvals.list
        printf '%s' '~{sub(sep("\n", betas), "'", "'\"'\"'")}' > betas.list
        printf '%s' '~{sub(sep("\n", beta_ses), "'", "'\"'\"'")}' > beta_ses.list
        printf '%s' '~{sub(sep("\n", afs), "'", "'\"'\"'")}' > afs.list
        python3 -u - '~{sub(prefix, "'", "'\"'\"'")}' ~{chromosome_count} ~{return_dense} <<'PY'
        from datetime import datetime, timezone
        from pathlib import Path
        import json
        import re
        import sys
        import pyarrow as pa
        import pyarrow.parquet as pq

        prefix, count, dense = sys.argv[1:]
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', prefix):
            sys.exit('Input error: prefix must be a filename prefix.')
        if int(count) < 1:
            sys.exit('Merge error: no chromosome tasks.')
        def log(message):
            stamp = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
            print(f'[{stamp}] stage=merge_trans {message}', file=sys.stderr)

        groups = [('pairs', 'pairs')] if dense == 'false' else [
            ('pvals', 'pval'), ('betas', 'beta'), ('beta_ses', 'beta_se'), ('afs', 'af')]
        for name, suffix in groups:
            paths = Path(name+'.list').read_text().splitlines()
            if len(paths) != int(count):
                sys.exit(f'Merge error: expected {count} {suffix} files; received {len(paths)}.')
            writer = None
            total = 0
            try:
                for path in paths:
                    if path.startswith(('gs://', 's3://', 'http://', 'https://')) or not Path(path).is_file():
                        sys.exit('Input localization error: unreadable result '+path)
                    parquet = pq.ParquetFile(path)
                    schema = parquet.schema_arrow
                    if dense == 'false':
                        pandas_metadata = json.loads((schema.metadata or {}).get(b'pandas', b'{}'))
                        indexes = [name for name in pandas_metadata.get('index_columns', []) if isinstance(name, str)]
                        schema = pa.schema([field for field in schema if field.name not in indexes])
                    if writer is None:
                        # Sparse row indexes restart on each chromosome. Dense files
                        # retain their physical variant_id index column and metadata.
                        output_schema = schema
                        writer = pq.ParquetWriter(prefix+'.trans_qtl_'+suffix+'.parquet', output_schema)
                    if not schema.equals(output_schema, check_metadata=False):
                        sys.exit('Merge error: incompatible result schema '+path)
                    for batch in parquet.iter_batches(batch_size=65536):
                        if dense == 'false':
                            batch = batch.select(output_schema.names)
                        table = pa.Table.from_batches([batch], schema=output_schema)
                        writer.write_table(table)
                        total += batch.num_rows
                    log(f'Merged {suffix}: {path}; rows={parquet.metadata.num_rows}.')
            finally:
                if writer is not None:
                    writer.close()
            log(f'Merge completed: {suffix}; rows={total}.')
        PY
    >>>
    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "~{auxiliary_memory}GB"
        cpu: 2
        disks: "local-disk ~{disk_space} HDD"
        preemptible: num_preempt
    }
    output {
        File? trans_qtl = prefix + ".trans_qtl_pairs.parquet"
        File? trans_qtls_pval = prefix + ".trans_qtl_pval.parquet"
        File? trans_qtl_beta = prefix + ".trans_qtl_beta.parquet"
        File? trans_qtl_beta_se = prefix + ".trans_qtl_beta_se.parquet"
        File? trans_qtl_af = prefix + ".trans_qtl_af.parquet"
    }
}

workflow tensorqtl_trans_workflow {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        File? interaction_file
        String prefix
        Float maf_threshold
        Float? fdr
        Boolean return_dense
        Float pval_threshold = 0.00001
        Int batch_size = 1000
        Int memory = 120
        Int disk_space
        Int num_threads = 32
        Int num_gpus = 1
        Int num_preempt
        Int split_memory = 8
        Int split_threads = 4
        Int auxiliary_memory = 16
    }
    call prepare_samples {
        input: phenotype_bed=phenotype_bed, covariates=covariates,
            interaction_file=interaction_file, prefix=prefix, return_dense=return_dense,
            disk_space=disk_space, auxiliary_memory=auxiliary_memory, num_preempt=num_preempt
    }
    call split_chromosomes {
        input: plink_pgen=plink_pgen, plink_pvar=plink_pvar, plink_psam=plink_psam,
            disk_space=disk_space, split_memory=split_memory, split_threads=split_threads,
            num_preempt=num_preempt
    }
    scatter (index in range(length(split_chromosomes.pgens))) {
        call tensorqtl_trans {
            input: plink_pgen=split_chromosomes.pgens[index],
                plink_pvar=split_chromosomes.pvars[index], plink_psam=split_chromosomes.psams[index],
                phenotype_bed=prepare_samples.bed, covariates=prepare_samples.cov,
                interaction_file=prepare_samples.interaction, prefix="~{prefix}.chr_~{index}",
                maf_threshold=maf_threshold, fdr=fdr, return_dense=return_dense,
                pval_threshold=pval_threshold, batch_size=batch_size, memory=memory,
                disk_space=disk_space, num_threads=num_threads, num_gpus=num_gpus,
                num_preempt=num_preempt
        }
    }
    call merge_trans {
        input: pairs=select_all(tensorqtl_trans.trans_qtl), pvals=select_all(tensorqtl_trans.trans_qtls_pval),
            betas=select_all(tensorqtl_trans.trans_qtl_beta), beta_ses=select_all(tensorqtl_trans.trans_qtl_beta_se),
            afs=select_all(tensorqtl_trans.trans_qtl_af), chromosome_count=length(split_chromosomes.pgens),
            prefix=prefix, return_dense=return_dense, disk_space=disk_space,
            auxiliary_memory=auxiliary_memory, num_preempt=num_preempt
    }
    output {
        File? trans_qtl = merge_trans.trans_qtl
        File? trans_qtls_pval = merge_trans.trans_qtls_pval
        File? trans_qtl_beta = merge_trans.trans_qtl_beta
        File? trans_qtl_beta_se = merge_trans.trans_qtl_beta_se
        File? trans_qtl_af = merge_trans.trans_qtl_af
        Array[String] chromosomes = split_chromosomes.chromosomes
    }
}
