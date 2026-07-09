"""Path bootstrap for standalone probes.

Adds the repository and scrutiny roots, imports the canonical
``interference.*`` modules, and maps historical short names to those objects.
"""
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))          # tests/probes
_tests = _os.path.dirname(_here)                              # tests/
_repo = _os.path.dirname(_tests)

for _p in (_os.path.join(_tests, "audit"),
           _os.path.join(_tests, "substrate"),
           _tests, _here, _repo):
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

from interference import (couplings, dark_sector, gravity, higgs_ew,
                          masses, mixing, reference, root, wiring, words)
from interference import _paths
from interference.searchlight import (embedding_uniqueness, octonions,
                                      protected_forgetting)

_compat_modules = {
    "_paths": _paths,
    "wiring": wiring,
    "root": root,
    "words": words,
    "masses": masses,
    "mixing": mixing,
    "couplings": couplings,
    "higgs_ew": higgs_ew,
    "gravity": gravity,
    "dark_sector": dark_sector,
    "reference": reference,
    "octonions": octonions,
    "embedding_uniqueness": embedding_uniqueness,
    "protected_forgetting": protected_forgetting,
}
for _name, _module in _compat_modules.items():
    _sys.modules[_name] = _module
