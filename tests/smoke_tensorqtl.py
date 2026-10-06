"""CPU smoke test of the rendered WDL command with the real tensorQTL engine."""
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd
import pgenlib
import WDL

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
        f'1\t{i+1}\tv{i}\tA\tG\n' for i in range(3)))
    (root/'input.psam').write_text('#IID\n' + '\n'.join(samples) + '\n')
    phenotype = pd.DataFrame([['2', 100, 101, 'GENE', *expression]],
                              columns=['#chr','start','end','phenotype_id',*samples])
    covariates = pd.DataFrame([rng.normal(size=n)], index=['PC1'], columns=samples)
    interactions = pd.Series(interaction, index=samples)
    phenotype.to_csv(root/'expression.bed', sep='\t', index=False)
    covariates.to_csv(root/'covariates.tsv', sep='\t')
    interactions.to_csv(root/'interaction.tsv', sep='\t', header=False)
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
    task = WDL.load(str(REPO/'tensorQTL_trans.wdl')).tasks[0]
    common = dict(plink_pgen=str(root/'input.pgen'), plink_pvar=str(root/'input.pvar'),
                  plink_psam=str(root/'input.psam'), phenotype_bed=str(root/'expression.bed'),
                  covariates=str(root/'covariates.tsv'), maf_threshold=0.05,
                  pval_threshold=1.0, memory=4, disk_space=10, num_threads=1, num_gpus=0, num_preempt=0)
    interaction_results = {}
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
        env = WDL.values_from_json(inputs, task.available_inputs)
        stdlib = WDL.StdLib.Base(task.effective_wdl_version)
        for decl in task.inputs:
            if decl.name not in [binding.name for binding in env]:
                env = env.bind(decl.name, decl.expr.eval(env, stdlib) if decl.expr else WDL.Value.Null())
        command = task.command.eval(env, stdlib).value
        subprocess.run(['bash','-c',command], cwd=work, check=True)
        if mode == 'dense':
            for suffix in ['pval','beta','beta_se','af']:
                assert (work/f'{mode}.trans_qtl_{suffix}.parquet').is_file()
        else:
            results = pd.read_parquet(work/f'{mode}.trans_qtl_pairs.parquet')
            assert set(results['variant_id']) == {'v0','v1','v2'}
            if mode in ['interaction', 'intersection', 'reference']:
                hit = results.loc[results['variant_id']=='v0'].iloc[0]
                assert hit['pval_gi'] < 1e-10, hit
                interaction_results[mode] = results
            if mode == 'intersection':
                aligned = pd.read_csv(work/'phenotype.aligned.bed.gz', sep='\t')
                assert list(aligned.columns[4:]) == retained
    pd.testing.assert_frame_equal(interaction_results['intersection'], interaction_results['reference'])
    print('Real tensorQTL smoke passed: interaction, three-file intersection versus reference, '
          'ordinary sparse, ordinary dense.')
