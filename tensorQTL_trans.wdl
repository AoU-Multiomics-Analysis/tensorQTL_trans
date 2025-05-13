task tensorqtl_trans {

    File plink_pgen
    File plink_pvar
    File plink_psam

    File phenotype_bed
    File covariates
    String prefix

    Float maf_threshold
    Float? fdr
    Boolean return_dense

    Int memory
    Int disk_space
    Int num_threads
    Int num_gpus
    Int num_preempt

    command {
        set -euo pipefail
        plink_base=$(echo "${plink_pgen}" | rev | cut -f 2- -d '.' | rev)
        python3 -m tensorqtl \
            $plink_base ${phenotype_bed} ${prefix} \
            --mode trans \
            --covariates ${covariates} \
            ${if return_dense then "--return_dense" else ""} \
            ${"--fdr " + fdr} \
            ${"--maf_threshold " + maf_threshold} 
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
        gpuType: "nvidia-tesla-p100"
        gpuCount: "${num_gpus}"
        zones: ["us-central1-c"]
    }

    output {
        File? trans_qtl_pairs = "${prefix}.trans_qtl_pairs.parquet"
        File? trans_qtls_pval = "${prefix}.trans_qtl_pval.parquet"
        File? trans_qtl_beta = "${prefix}.trans_qtl_beta.parquet"
        File? trans_qtl_beta_se = "${prefix}.trans_qtl_beta_se.parquet"
        File? trans_qtl_af = "${prefix}.trans_qtl_af.parquet"
    }    
    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_trans_workflow {
    call tensorqtl_trans
}
