"""The reusable optical viewer does not require the QE host application."""

import os
import subprocess
import sys


def test_spectrum_widgets_import_without_qe(tmp_path):
    code = """
import importlib.abc
import sys

class NoQe(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "aiidalab_qe" or fullname.startswith("aiidalab_qe."):
            raise ModuleNotFoundError("The standalone viewer must not import QE")

sys.meta_path.insert(0, NoQe())
from aiidalab_qe_vibroscopy.app.widgets.ramanmodel import RamanModel
from aiidalab_qe_vibroscopy.app.widgets.ramanwidget import RamanWidget
model = RamanModel(spectrum_type="Raman")
assert model.spectrum_type == "Raman"
assert not any(name.startswith("aiidalab_qe.") for name in sys.modules)
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=os.environ | {"AIIDA_PATH": str(tmp_path / "state")},
    )
    assert result.returncode == 0, result.stdout + result.stderr
