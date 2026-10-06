"""Render task commands with typed WDL inputs, as used by command smoke tests."""
from pathlib import Path
import tempfile
import WDL

REPO = Path(__file__).resolve().parents[1]
_COMMAND_FILES = tempfile.TemporaryDirectory(prefix='wdl-smoke-command-files-')

class CommandStdLib(WDL.StdLib.Base):
    """Keep command-rendered files readable until test subprocesses finish."""
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        return filename

def render_task(name, inputs, local_paths=None):
    task = next(t for t in WDL.load(str(REPO/'tensorQTL_trans.wdl')).tasks if t.name == name)
    values = WDL.values_from_json({k: v for k, v in inputs.items()
                                  if k in {b.name for b in task.available_inputs}}, task.available_inputs)
    if local_paths:
        values = WDL.Value.rewrite_env_paths(values, lambda value: local_paths[value.value])
    stdlib = CommandStdLib(task.effective_wdl_version, write_dir=_COMMAND_FILES.name)
    for decl in task.inputs:
        if decl.name not in {b.name for b in values}:
            values = values.bind(decl.name, decl.expr.eval(values, stdlib) if decl.expr else WDL.Value.Null())
    return task.command.eval(values, stdlib).value
