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
samples = [f'S{i}' for i in range(n)]
genotypes = rng.binomial(2, 0.35, size=(3, n)).astype(np.int8)
interaction = rng.uniform(0.05, 0.5, n)
expression = 2 + 8*genotypes[0]*interaction + rng.normal(0, 0.1, n)

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    with pgenlib.PgenWriter(os.fsencode(root/'input.pgen'), sample_ct=n, variant_ct=3) as writer:
        for row in genotypes:
            writer.append(row)
    (root/'input.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'1\t{i+1}\tv{i}\tA\tG\n' for i in range(3)))
    (root/'input.psam').write_text('#IID\n' + '\n'.join(samples) + '\n')
    pd.DataFrame([['2', 100, 101, 'GENE', *expression]],
                 columns=['#chr','start','end','phenotype_id',*samples]).to_csv(root/'expression.bed', sep='\t', index=False)
    pd.DataFrame([rng.normal(size=n)], index=['PC1'], columns=samples).to_csv(root/'covariates.tsv', sep='\t')
    pd.Series(interaction, index=samples).to_csv(root/'interaction.tsv', sep='\t', header=False)
    task = WDL.load(str(REPO/'tensorQTL_trans.wdl')).tasks[0]
    common = dict(plink_pgen=str(root/'input.pgen'), plink_pvar=str(root/'input.pvar'),
                  plink_psam=str(root/'input.psam'), phenotype_bed=str(root/'expression.bed'),
                  covariates=str(root/'covariates.tsv'), maf_threshold=0.05,
                  pval_threshold=1.0, memory=4, disk_space=10, num_threads=1, num_gpus=0, num_preempt=0)
    for mode in ['interaction', 'ordinary', 'dense']:
        work = root/mode
        work.mkdir()
        inputs = dict(common, prefix=mode, return_dense=mode=='dense')
        if mode == 'interaction':
            inputs['interaction_file'] = str(root/'interaction.tsv')
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
            if mode == 'interaction':
                hit = results.loc[results['variant_id']=='v0'].iloc[0]
                assert hit['pval_gi'] < 1e-10, hit
    print('Real tensorQTL smoke passed: interaction, ordinary sparse, ordinary dense.')
