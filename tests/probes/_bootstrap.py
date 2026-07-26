"""Path bootstrap for the probes.

Puts interference/ (flat: the derivation) and its searchlight/
subdirectory on sys.path via the theory's own _paths, then adds the
scrutiny directories: tests/substrate/ (the NLS conjecture) and
tests/audit/ (the committed grammar and the closure search).  A
probe's bare-name imports (root, substrate3d, edge_grammar, ...)
resolve either way, both under pytest and standalone
(python3 <probe>.py).
"""
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))          # tests/probes
_tests = _os.path.dirname(_here)                              # tests/
_repo = _os.path.dirname(_tests)
_iface = _os.path.join(_repo, "interference")

for _p in (_os.path.join(_tests, "audit"),
           _os.path.join(_tests, "substrate"),
           _tests, _here, _iface):
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import _paths  # noqa: E402,F401  (adds interference/ and searchlight/)
