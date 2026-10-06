"""CPU smoke test of the rendered WDL command with the real tensorQTL engine."""
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd
import pgenlib
from task_runner import render_task

REPO = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(2026)
n = 120
# Numeric IDs with leading zeros reproduce tensorQTL's interaction-index inference bug.
samples = [f'{1000+i:08d}' for i in range(n)]
genotypes = rng.binomial(2, 0.35, size=(3, n)).astype(np.int8)
interaction = rng.uniform(0.05, 0.5, n)
expression = 2 + 8*genotypes[0]*interaction + rng.normal(0, 0.1, n)

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    with pgenlib.PgenWriter(os.fsencode(root/'input.pgen'), sample_ct=n, variant_ct=3) as writer:
        for row in genotypes:
            writer.append_biallelic(row)
    (root/'input.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'{chrom}\t{i+1}\tv{i}\tA\tG\n' for i, chrom in enumerate(['chr1', 'chr2', 'chrX'])))
    (root/'input.psam').write_text('#IID\n' + '\n'.join(samples) + '\n')
    phenotype = pd.DataFrame([['chr9', 100, 101, 'GENE', *expression],
                              ['chr1', 1, 2, 'CIS_GENE', *(3+genotypes[1]*interaction+rng.normal(0,0.2,n))]],
                              columns=['#chr','start','end','phenotype_id',*samples])
    covariates = pd.DataFrame([rng.normal(size=n)], index=['PC1'], columns=samples)
    interactions = pd.Series(interaction, index=samples)
    phenotype.to_csv(root/'expression.bed', sep='\t', index=False)
    covariates.to_csv(root/'covariates.tsv', sep='\t')
    interactions.to_csv(root/'interaction.tsv', sep='\t', header=False)
    # Keep ordinary modes on their direct CLI path, using IDs without leading zeros.
    ordinary_samples = [str(int(sample)) for sample in samples]
    rename_samples = dict(zip(samples, ordinary_samples))
    phenotype.rename(columns=rename_samples).to_csv(root/'ordinary.bed', sep='\t', index=False)
    covariates.rename(columns=rename_samples).to_csv(root/'ordinary.covariates.tsv', sep='\t')
    (root/'ordinary.psam').write_text('#IID\n'+'\n'.join(ordinary_samples)+'\n')
    # The three sample sets differ, and both TSV files use a different order.
    phenotype.to_csv(root/'partial.bed.gz', sep='\t', index=False)
    partial_covariates = covariates.loc[:, samples[10:][::-1]].copy()
    partial_covariates['covariate_only'] = 0.0
    partial_covariates.to_csv(root/'partial.covariates.tsv', sep='\t')
    partial_interactions = interactions.loc[samples[:110][::-1]].copy()
    partial_interactions.loc['interaction_only'] = 0.25
    partial_interactions.to_csv(root/'partial.interaction.tsv', sep='\t', header=False)
    retained = samples[10:110]
    phenotype.loc[:, ['#chr', 'start', 'end', 'phenotype_id', *retained]].to_csv(
        root/'reference.bed', sep='\t', index=False)
    covariates.loc[:, retained].to_csv(root/'reference.covariates.tsv', sep='\t')
    interactions.loc[retained].to_csv(root/'reference.interaction.tsv', sep='\t', header=False)
    common = dict(plink_pgen=str(root/'input.pgen'), plink_pvar=str(root/'input.pvar'),
                  plink_psam=str(root/'input.psam'), phenotype_bed=str(root/'expression.bed'),
                  covariates=str(root/'covariates.tsv'), maf_threshold=0.05,
                  pval_threshold=1.0, memory=4, disk_space=10, num_threads=1, num_gpus=0, num_preempt=0)
    interaction_results = {}
    mode_inputs = {}
    mode_results = {}
    for mode in ['interaction', 'intersection', 'reference', 'ordinary', 'dense']:
        work = root/mode
        work.mkdir()
        inputs = dict(common, prefix=mode, return_dense=mode=='dense')
        if mode == 'interaction':
            inputs['interaction_file'] = str(root/'interaction.tsv')
        elif mode in ['intersection', 'reference']:
            stem = 'partial' if mode == 'intersection' else 'reference'
            inputs.update(phenotype_bed=str(root/(stem+'.bed.gz' if mode == 'intersection' else stem+'.bed')),
                          covariates=str(root/(stem+'.covariates.tsv')),
                          interaction_file=str(root/(stem+'.interaction.tsv')))
        else:
            inputs.update(plink_psam=str(root/'ordinary.psam'),
                          phenotype_bed=str(root/'ordinary.bed'),
                          covariates=str(root/'ordinary.covariates.tsv'))
        subprocess.run(['bash','-c',render_task('prepare_samples',inputs)],cwd=work,check=True)
        if 'interaction_file' in inputs:
            inputs.update(phenotype_bed=str(work/'phenotype.aligned.bed.gz'),
                          covariates=str(work/'covariates.aligned.tsv'),
                          interaction_file=str(work/'prepared.interaction.tsv'))
        mode_inputs[mode] = inputs.copy()
        subprocess.run(['bash','-c',render_task('tensorqtl_trans',inputs)],cwd=work,check=True)
        if mode == 'dense':
            for suffix in ['pval','beta','beta_se','af']:
                assert (work/f'{mode}.trans_qtl_{suffix}.parquet').is_file()
        else:
            results = pd.read_parquet(work/f'{mode}.trans_qtl_pairs.parquet')
            assert set(results['variant_id']) == {'v0','v1','v2'}
            assert not ((results.variant_id=='v0') & (results.phenotype_id=='CIS_GENE')).any()
            mode_results[mode] = results
            if mode in ['interaction', 'intersection', 'reference']:
                hit = results.loc[(results['variant_id']=='v0') & (results['phenotype_id']=='GENE')].iloc[0]
                assert hit['pval_gi'] < 1e-10, hit
                interaction_results[mode] = results
            if mode == 'intersection':
                aligned = pd.read_csv(work/'phenotype.aligned.bed.gz', sep='\t')
                assert list(aligned.columns[4:]) == retained
    pd.testing.assert_frame_equal(interaction_results['intersection'], interaction_results['reference'])
    # Execute the complete task chain with real PLINK2 and tensorQTL, without
    # submitting a Terra job. Use the exact same prepared sample files in each run.
    split = root/'split'
    split.mkdir()
    subprocess.run(['bash','-c',render_task('split_chromosomes',common)],cwd=split,check=True)
    chromosome_files = sorted((split/'chromosomes').glob('*.pgen'))
    assert len(chromosome_files) == 3
    assert (split/'chromosomes.txt').read_text().splitlines() == ['chr1','chr2','chrX']
    for i,path in enumerate(chromosome_files):
        with pgenlib.PgenReader(os.fsencode(path)) as reader:
            observed = np.empty(n,dtype=np.int8)
            reader.read(0,observed)
            np.testing.assert_array_equal(observed,genotypes[i])
        assert path.with_suffix('.psam').read_text() == (root/'input.psam').read_text()
    for mode in ['interaction','intersection','ordinary','dense','mixed_hits','no_hits']:
        inputs = mode_inputs['intersection' if mode in ['mixed_hits','no_hits'] else mode].copy()
        dense = mode=='dense'
        if mode in ['mixed_hits','no_hits']:
            inputs['pval_threshold'] = 1e-30 if mode=='mixed_hits' else 1e-100
        groups = {key:[] for key in ['pairs','pvals','betas','beta_ses','afs']}
        for i,path in enumerate(chromosome_files):
            work = root/f'scatter_{mode}_{i}'
            work.mkdir()
            shard_inputs = dict(inputs,plink_pgen=str(path),plink_pvar=str(path.with_suffix('.pvar')),
                                plink_psam=str(path.with_suffix('.psam')),prefix='shard')
            if mode in ['ordinary','dense']:
                shard_inputs['plink_psam'] = str(root/'ordinary.psam')
            subprocess.run(['bash','-c',render_task('tensorqtl_trans',shard_inputs)],cwd=work,check=True)
            for key,suffix in [('pairs','pairs'),('pvals','pval'),('betas','beta'),('beta_ses','beta_se'),('afs','af')]:
                output = work/f'shard.trans_qtl_{suffix}.parquet'
                if output.is_file():
                    groups[key].append(str(output))
        work = root/f'merged_{mode}'
        work.mkdir()
        merge_inputs = dict(common,**groups,chromosome_count=3,prefix='merged',return_dense=dense)
        subprocess.run(['bash','-c',render_task('merge_trans',merge_inputs)],cwd=work,check=True)
        if dense:
            for suffix in ['pval','beta','beta_se','af']:
                expected = pd.read_parquet(root/f'dense/dense.trans_qtl_{suffix}.parquet')
                observed = pd.read_parquet(work/f'merged.trans_qtl_{suffix}.parquet')
                pd.testing.assert_frame_equal(observed.sort_index(),expected.sort_index(),rtol=1e-3,atol=1e-6)
        else:
            observed = pd.read_parquet(work/'merged.trans_qtl_pairs.parquet')
            if mode=='no_hits':
                assert observed.empty
            elif mode=='mixed_hits':
                assert set(observed.variant_id) == {'v0'}
                assert observed.pval_gi.max() < 1e-30
            else:
                def order(frame):
                    return frame.sort_values(['variant_id','phenotype_id']).reset_index(drop=True)
                pd.testing.assert_frame_equal(order(observed),order(mode_results[mode]),rtol=1e-3,atol=1e-6)
    print('Real PLINK2/tensorQTL smoke passed: chromosome genotype integrity, '
          'merged versus whole-genome interaction, sample intersection, ordinary sparse, '
          'dense indexes, cis exclusion, mixed-hit and no-hit interaction results.')
