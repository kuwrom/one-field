"""Regression tests for the derivation's single canonical module graph."""

import os
import subprocess
import sys


def test_package_imports_do_not_create_bare_module_duplicates():
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    code = """
import sys
from interference import couplings, dark_sector, gravity, masses, mixing, root
from interference.searchlight import octonions

for name in ("wiring", "root", "masses", "mixing", "couplings",
             "higgs_ew", "gravity", "dark_sector", "octonions"):
    assert name not in sys.modules, name
assert masses.derive.__module__ == "interference.masses"
assert root.WEB.__class__.__module__ == "interference.root"
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=repo,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
