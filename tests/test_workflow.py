"""Render and execute the real task command without a GPU or Docker daemon."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import WDL

ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def run_task(self, interaction=None, dense=False, cloud=False, covariate='PC1\t0\t1\t2\n'):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {}
            # Separate localization directories, with spaces and shell metacharacters.
            for key, name in [('plink_pgen', 'geno.pgen'), ('plink_pvar', 'different.pvar'),
                              ('plink_psam', 'other.psam'), ('phenotype_bed', 'expression.bed'),
                              ('covariates', 'covariates.tsv')]:
                folder = root / (key + " space '$literal")
                folder.mkdir()
                path = folder / name
                path.write_text('placeholder\n')
                paths[key] = str(path)
            Path(paths['phenotype_bed']).write_text('#chr\tstart\tend\tphenotype_id\tS1\tS2\tS3\n9\t0\t1\tGENE\t1\t2\t3\n')
            Path(paths['covariates']).write_text('ID\tS1\tS2\tS3\n' + covariate)
            if interaction is not None:
                ipath = root / "interaction space '$literal.tsv"
                ipath.write_text(interaction)
                paths['interaction_file'] = 'gs://bucket/reference.tsv' if cloud else str(ipath)
            inputs = dict(paths, prefix='test', maf_threshold=0.05, return_dense=dense,
                          memory=4, disk_space=10, num_threads=1, num_gpus=0, num_preempt=0)
            task = WDL.load(str(ROOT / 'tensorQTL_trans.wdl')).tasks[0]
            self.assertIn('interaction_file', [b.name for b in task.available_inputs],
                          'Task must expose the optional interaction File input')
            env = WDL.values_from_json(inputs, task.available_inputs)
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
                'import sys,json,pathlib\n'
                'pathlib.Path("argv.json").write_text(json.dumps(sys.argv[1:]))\n'
                'assert pathlib.Path(sys.argv[1]+".pvar").read_text()=="placeholder\\n"\n'
                'assert pathlib.Path(sys.argv[1]+".psam").read_text()=="placeholder\\n"\n'
                'if "--interaction" in sys.argv:\n'
                ' p=pathlib.Path(sys.argv[sys.argv.index("--interaction")+1])\n'
                ' assert p.read_text()=="S1\\t0.1\\nS2\\t0.3\\nS3\\t0.2\\n"\n')
            result = subprocess.run(['bash', '-c', command], cwd=root, text=True,
                                    capture_output=True, env={**os.environ, 'PYTHONPATH':str(root)})
            args = json.loads((root/'argv.json').read_text()) if (root/'argv.json').exists() else None
            return result, args

    def test_optional_input_absent(self):
        result, args = self.run_task()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('--interaction', args)

    def test_localized_interaction_is_aligned_and_passed(self):
        result, args = self.run_task('S3\t0.2\nS1\t0.1\nS2\t0.3\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--interaction', args)

    def test_reject_invalid_interactions_before_engine(self):
        for data in ['S1\t0.1\nS2\t0.3\n', 'S1\t0.1\nS1\t0.3\nS3\t0.2\n',
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
        result, args = self.run_task('S1\t0.1\nS2\t0.3\nS3\t0.2\n', cloud=True)
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
