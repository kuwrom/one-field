"""Repository-root path bootstrap for standalone scrutiny scripts.

Canonical derivation modules use the ``interference.*`` namespace.
"""
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))
_repo = _os.path.dirname(_here)
if _repo not in _sys.path:
    _sys.path.insert(0, _repo)
