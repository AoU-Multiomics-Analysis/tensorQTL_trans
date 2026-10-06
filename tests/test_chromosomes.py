import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import pandas as pd
from task_runner import render_task

class ChromosomeTests(unittest.TestCase):
    def run_command(self, command, root):
        return subprocess.run(['bash', '-c', command], cwd=root, text=True, capture_output=True,
                              env={**os.environ, 'PATH':str(root)+os.pathsep+
                                   str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH']})

    def test_split_preserves_all_chromosomes_and_metadata_with_localized_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            local = {}
            inputs = dict(disk_space=10, split_memory=4, split_threads=1, num_preempt=0)
            contents = {
                'plink_pgen': 'fake binary',
                'plink_pvar': '##source=test\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t1\tv1\tA\tG\nchr2\t2\tv2\tC\tT\nchrX\t3\tvX\tA\tT\n',
                'plink_psam': '#IID\n001\n002\n'}
            for key, content in contents.items():
                path = root/(key+" space '$literal")
                path.write_text(content)
                uri = 'gs://bucket/'+key
                inputs[key] = uri
                local[uri] = str(path)
            # PLINK itself is covered by the real GitHub Actions smoke. This double
            # deliberately normalizes chromosome labels, as real PLINK does.
            (root/'plink2').write_text('#!'+sys.executable+'\n'+'''import pathlib,sys
args=sys.argv[1:]
def value(flag): return args[args.index(flag)+1]
out=value('--out')
chrom=value('--chr')
rows=pathlib.Path(value('--pvar')).read_text().splitlines()
selected=[row for row in rows if not row.startswith('#') and row.split()[0]==chrom]
pathlib.Path(out+'.pgen').write_text('fake binary')
pathlib.Path(out+'.pvar').write_text('#CHROM\\tPOS\\tID\\tREF\\tALT\\n'+'\\n'.join(row.replace('chr','',1) for row in selected)+'\\n')
pathlib.Path(out+'.psam').write_text(pathlib.Path(value('--psam')).read_text())
''')
            (root/'plink2').chmod(0o755)
            result = self.run_command(render_task('split_chromosomes', inputs, local), root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root/'chromosomes.txt').read_text().splitlines(), ['chr1','chr2','chrX'])
            for i, chrom in enumerate(['chr1','chr2','chrX'], 1):
                stem = root/f'chromosomes/chr_{i:06d}'
                self.assertEqual(Path(str(stem)+'.psam').read_text(), contents['plink_psam'])
                self.assertEqual(Path(str(stem)+'.pvar').read_text().splitlines()[2].split()[0], chrom)
            self.assertIn('Split completed: 3 chromosomes', result.stderr)

    def test_split_rejects_unlocalized_files_and_empty_pvar(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dummy = root/'dummy'
            dummy.write_text('#CHROM\tPOS\tID\tREF\tALT\n')
            common = dict(plink_pgen=str(dummy), plink_pvar=str(dummy), plink_psam=str(dummy),
                          disk_space=10, num_preempt=0)
            for key in ['plink_pgen','plink_pvar','plink_psam']:
                result = self.run_command(render_task('split_chromosomes', {**common, key:'gs://bucket/input'}),root)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('localization error',result.stderr)
            result = self.run_command(render_task('split_chromosomes',common),root)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('no variants',result.stderr)

    def merge(self, root, frames, dense=False, localize=False):
        paths = []
        mapping = {}
        for i, frame in enumerate(frames):
            path = root/f"input {i} '$literal.parquet"
            frame.to_parquet(path)
            uri = f'gs://bucket/{i}.parquet'
            paths.append(uri if localize else str(path))
            mapping[uri] = str(path)
        inputs = dict(pairs=[] if dense else paths, pvals=paths if dense else [],
                      betas=paths if dense else [], beta_ses=paths if dense else [], afs=paths if dense else [],
                      chromosome_count=len(frames), prefix='merged', return_dense=dense,
                      disk_space=10, num_preempt=0)
        return self.run_command(render_task('merge_trans',inputs,mapping if localize else None),root)

    def test_sparse_merge_handles_empty_chromosomes_and_localization(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Cis exclusion leaves a physical pandas row-index column in some files.
            first = pd.DataFrame({'variant_id':['v1'], 'phenotype_id':['gene'], 'pval_gi':[0.001]},index=[7])
            empty = pd.DataFrame({'variant_id':pd.Series(dtype='string'),
                                  'phenotype_id':pd.Series(dtype='string'), 'pval_gi':pd.Series(dtype='float64')})
            last = first.assign(variant_id='v2')
            result = self.merge(root,[empty,first,empty,last],localize=True)
            self.assertEqual(result.returncode,0,result.stderr)
            merged = pd.read_parquet(root/'merged.trans_qtl_pairs.parquet')
            pd.testing.assert_frame_equal(merged,pd.concat([first,last],ignore_index=True))

    def test_all_empty_sparse_results_are_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty = pd.DataFrame({'variant_id':pd.Series(dtype='string'), 'pval_gi':pd.Series(dtype='float64')})
            result = self.merge(root,[empty,empty])
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(len(pd.read_parquet(root/'merged.trans_qtl_pairs.parquet')),0)

    def test_dense_merge_preserves_variant_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frames = [pd.DataFrame({'gene':[0.01]},index=pd.Index([v],name='variant_id')) for v in ['v1','v2']]
            result = self.merge(root,frames,dense=True,localize=True)
            self.assertEqual(result.returncode,0,result.stderr)
            for suffix in ['pval','beta','beta_se','af']:
                pd.testing.assert_frame_equal(pd.read_parquet(root/f'merged.trans_qtl_{suffix}.parquet'),pd.concat(frames))

    def test_merge_rejects_schema_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.merge(Path(tmp),[pd.DataFrame({'variant_id':['v1']}), pd.DataFrame({'wrong':['v2']})])
            self.assertNotEqual(result.returncode,0)
            self.assertIn('incompatible result schema',result.stderr)

    def test_merge_rejects_missing_chromosome_outputs_and_cloud_uris(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            common=dict(pairs=[],pvals=[],betas=[],beta_ses=[],afs=[],chromosome_count=2,
                        prefix='merged',return_dense=False,disk_space=10,num_preempt=0)
            result=self.run_command(render_task('merge_trans',common),root)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('expected 2 pairs files',result.stderr)
            common.update(pairs=['gs://bucket/a','gs://bucket/b'])
            result=self.run_command(render_task('merge_trans',common),root)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('localization error',result.stderr)
