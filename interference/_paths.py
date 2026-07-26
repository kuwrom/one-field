"""Path bootstrap for the theory's bare-name imports.

The repository is read-and-run, not a package: modules import each
other by bare name (root, masses, wiring, ...).  The derivation is
flat, so interference/ itself carries it; searchlight/ is the one
subdirectory, set aside because it is logically dispensable (the
octonions located the coordinates; the engine does not import them).
Every entry point inserts interference/ and then imports this
module.
"""
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))
for _d in ("searchlight", None):
    _p = _here if _d is None else _os.path.join(_here, _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
