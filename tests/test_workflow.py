"""Render and execute the real task command without a GPU or Docker daemon."""
import json
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import WDL

ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def run_task(self, interaction=None, dense=False, cloud=False, covariate='PC1\t0\t1\t2\n',
                 covariates_data=None, phenotype_data=None, compressed=False, localize=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {}
            # Separate localization directories, with spaces and shell metacharacters.
            for key, name in [('plink_pgen', 'geno.pgen'), ('plink_pvar', 'different.pvar'),
                              ('plink_psam', 'other.psam'),
                              ('phenotype_bed', 'expression.bed.gz' if compressed else 'expression.bed'),
                              ('covariates', 'covariates.tsv')]:
                folder = root / (key + " space '$literal")
                folder.mkdir()
                path = folder / name
                path.write_text('placeholder\n')
                paths[key] = str(path)
            phenotype_data = phenotype_data if phenotype_data is not None else (
                '#chr\tstart\tend\tphenotype_id\tS1\tS2\tS3\n9\t0\t1\tGENE\t1\t2\t3\n')
            opener = gzip.open if compressed else open
            with opener(paths['phenotype_bed'], 'wt') as handle:
                handle.write(phenotype_data)
            Path(paths['covariates']).write_text(
                covariates_data if covariates_data is not None else 'ID\tS1\tS2\tS3\n' + covariate)
            if interaction is not None:
                ipath = root / "interaction space '$literal.tsv"
                ipath.write_text(interaction)
                paths['interaction_file'] = str(ipath)
            if cloud:
                paths['interaction_file' if cloud is True else cloud] = 'gs://bucket/reference.tsv'
            inputs = dict(paths, prefix='test', maf_threshold=0.05, return_dense=dense,
                          memory=4, disk_space=10, num_threads=1, num_gpus=0, num_preempt=0)
            task = WDL.load(str(ROOT / 'tensorQTL_trans.wdl')).tasks[0]
            self.assertIn('interaction_file', [b.name for b in task.available_inputs],
                          'Task must expose the optional interaction File input')
            if localize:
                # Simulate Cromwell replacing typed cloud File inputs at command rendering.
                local_paths = {f'gs://bucket/{key}/{Path(path).name}': path
                               for key, path in paths.items()}
                inputs.update({key: f'gs://bucket/{key}/{Path(path).name}'
                               for key, path in paths.items()})
            env = WDL.values_from_json(inputs, task.available_inputs)
            if localize:
                env = WDL.Value.rewrite_env_paths(env, lambda value: local_paths[value.value])
            stdlib = WDL.StdLib.Base(task.effective_wdl_version)
            for decl in task.inputs:
                if decl.name not in [binding.name for binding in env]:
                    value = decl.expr.eval(env, stdlib) if decl.expr else WDL.Value.Null()
                    env = env.bind(decl.name, value)
            command = task.command.eval(env, stdlib).value
            # Replace only the external tensorQTL engine, retaining shell/localization/validation.
            package = root / 'tensorqtl'
            package.mkdir()
            (package / '__init__.py').write_text('')
            (package / '__main__.py').write_text(
                'import sys,json,pathlib,csv,gzip\n'
                'pathlib.Path("argv.json").write_text(json.dumps(sys.argv[1:]))\n'
                'assert pathlib.Path(sys.argv[1]+".pvar").read_text()=="placeholder\\n"\n'
                'assert pathlib.Path(sys.argv[1]+".psam").read_text()=="placeholder\\n"\n'
                'data={}\n'
                'for key,p in [("phenotype",sys.argv[2]),'
                '("covariates",sys.argv[sys.argv.index("--covariates")+1])]:\n'
                ' opener=gzip.open if p.endswith(".gz") else open\n'
                ' with opener(p,"rt") as f: data[key]=list(csv.reader(f,delimiter="\\t"))\n'
                'if "--interaction" in sys.argv:\n'
                ' with open(sys.argv[sys.argv.index("--interaction")+1]) as f:\n'
                '  data["interaction"]=list(csv.reader(f,delimiter="\\t"))\n'
                'pathlib.Path("mapping_inputs.json").write_text(json.dumps(data))\n')
            result = subprocess.run(['bash', '-c', command], cwd=root, text=True,
                                    capture_output=True, env={**os.environ, 'PYTHONPATH':str(root)})
            args = json.loads((root/'argv.json').read_text()) if (root/'argv.json').exists() else None
            result.mapping_inputs = (json.loads((root/'mapping_inputs.json').read_text())
                                     if (root/'mapping_inputs.json').exists() else None)
            return result, args

    def test_optional_input_absent(self):
        result, args = self.run_task()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('--interaction', args)
        self.assertEqual(result.mapping_inputs['phenotype'][0][4:], ['S1', 'S2', 'S3'])
        self.assertEqual(result.mapping_inputs['covariates'],
                         [['ID', 'S1', 'S2', 'S3'], ['PC1', '0', '1', '2']])

    def test_typed_cloud_files_are_localized_before_intersection(self):
        result, args = self.run_task('S2\t0.3\nS1\t0.1\n',
                                     covariate='PC1\t7\t7\t2\n', compressed=True, localize=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.mapping_inputs['phenotype'][0][4:], ['S1', 'S2'])
        self.assertEqual(result.mapping_inputs['covariates'],
                         [['ID', 'S1', 'S2'], ['PC1', '7', '7']])
        self.assertEqual(result.mapping_inputs['interaction'], [['S1', '0.1'], ['S2', '0.3']])

    def test_localized_interaction_is_aligned_and_passed(self):
        result, args = self.run_task('S3\t0.2\nS1\t0.1\nS2\t0.3\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--interaction', args)
        self.assertEqual(result.mapping_inputs['interaction'],
                         [['S1', '0.1'], ['S2', '0.3'], ['S3', '0.2']])

    def test_three_way_intersection_keeps_bed_order_and_values(self):
        phenotype = ('#chr\tstart\tend\tphenotype_id\tS3\tS1\tS2\tS4\tS7\n'
                     '9\t0\t1\tGENE\t30\t10\t20\t40\t70\n'
                     '9\t1\t2\tGENE2\t300\t100\t200\t400\t700\n')
        covariates = 'ID\tS5\tS2\tS3\tS1\nPC1\t50\t20\t30\t10\n'
        for compressed in [False, True]:
            with self.subTest(compressed=compressed):
                result, args = self.run_task('S4\t0.4\nS2\t0.3\nS1\t0.1\nS6\t0.6\nS3\t0.2\n',
                                             covariates_data=covariates,
                                             phenotype_data=phenotype, compressed=compressed)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.mapping_inputs['phenotype'], [
                    ['#chr', 'start', 'end', 'phenotype_id', 'S3', 'S1', 'S2'],
                    ['9', '0', '1', 'GENE', '30', '10', '20'],
                    ['9', '1', '2', 'GENE2', '300', '100', '200']])
                self.assertEqual(result.mapping_inputs['covariates'],
                                 [['ID', 'S3', 'S1', 'S2'], ['PC1', '30', '10', '20']])
                self.assertEqual(result.mapping_inputs['interaction'],
                                 [['S3', '0.2'], ['S1', '0.1'], ['S2', '0.3']])
                self.assertIn('retained=3', result.stderr)

    def test_missing_interaction_samples_are_removed(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.3\n', covariate='PC1\t7\t7\t2\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.mapping_inputs['phenotype'][0][4:], ['S1', 'S2'])

    def test_covariates_are_reordered_with_no_sample_loss(self):
        result, args = self.run_task('S3\t0.2\nS1\t0.1\nS2\t0.3\n',
                                     covariates_data='ID\tS3\tS1\tS2\nPC1\t2\t0\t1\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.mapping_inputs['covariates'],
                         [['ID', 'S1', 'S2', 'S3'], ['PC1', '0', '1', '2']])

    def test_reject_empty_intersection_before_engine(self):
        result, args = self.run_task('S4\t0.1\nS5\t0.3\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('no shared samples', result.stderr.lower())
        self.assertIsNone(args)

    def test_reject_interaction_constant_after_intersection(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.1\nS4\t0.3\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('vary', result.stderr)
        self.assertIsNone(args)

    def test_reject_duplicate_main_effect_after_intersection(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n',
                                     covariates_data='ID\tS3\tS1\tS2\tS4\nCD4\t20\t10\t30\t99\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('main effect', result.stderr)
        self.assertIsNone(args)

    def test_reject_ambiguous_covariate_samples(self):
        for covariates in ['ID\tS1\tS1\tS3\nPC1\t0\t1\t2\n',
                           'ID\tS1\t\tS3\nPC1\t0\t1\t2\n', '']:
            with self.subTest(covariates=covariates):
                result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n',
                                             covariates_data=covariates)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('covariate', result.stderr)
                self.assertIsNone(args)

    def test_reject_malformed_bed_rows_before_engine(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.3\n', phenotype_data=(
            '#chr\tstart\tend\tphenotype_id\tS1\tS2\tS3\n9\t0\t1\tGENE\t1\t2\n'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('BED rows', result.stderr)
        self.assertIsNone(args)

    def test_reject_invalid_interactions_before_engine(self):
        for data in ['S1\t0.1\nS1\t0.3\nS3\t0.2\n',
                     'S1\tNaN\nS2\t0.3\nS3\t0.2\n', 'S1\t1\t2\nS2\t2\t3\n',
                     'S1\t0.2\nS2\t0.2\nS3\t0.2\n', 'sample\tCD4\nS1\t0.1\nS2\t0.3\nS3\t0.2\n']:
            with self.subTest(data=data):
                result, args = self.run_task(data)
                self.assertNotEqual(result.returncode, 0)
                self.assertIsNone(args)

    def test_reject_duplicate_interaction_main_effect(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n',
                                     covariate='CD4\t0.1\t0.3\t0.2\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIsNone(args)

    def test_reject_unlocalized_cloud_input(self):
        for key in ['interaction_file', 'phenotype_bed', 'covariates']:
            with self.subTest(key=key):
                result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n', cloud=key)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('localization', result.stderr.lower())
                self.assertIsNone(args)

    def test_reject_dense_interaction_before_engine(self):
        result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n', dense=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIsNone(args)

    def test_dense_without_interaction(self):
        result, args = self.run_task(dense=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--return_dense', args)
        self.assertNotIn('--interaction', args)

    def test_no_workflow_scope_file_writes(self):
        doc = WDL.load(str(ROOT/'tensorQTL_trans.wdl'))
        forbidden = {'write_lines', 'write_tsv', 'write_map', 'write_json'}
        def visit(node):
            if isinstance(node, WDL.Expr.Apply):
                self.assertNotIn(node.function_name, forbidden)
            for child in node.children:
                visit(child)
        visit(doc.workflow)

if __name__ == '__main__':
    unittest.main()
