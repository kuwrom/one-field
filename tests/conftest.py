"""Shared fixtures for the scorecard test suite.

The full derivation chain is run once per test session, silently.  Every
test then asserts against the cached results, so the suite is fast and
the detailed console output never pollutes pytest.

Run with:
    pytest -q
"""

import contextlib
import io
import os
import sys

import pytest

# Add the repository root and tests/ so the canonical ``interference.*``
# package and the scrutiny helpers resolve consistently.
_repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)

# Add tests/ itself so the mechanism probes import as `probes.*`, plus
# the two scrutiny directories: substrate/ (the NLS substrate) and
# audit/ (grammar + closure search), whose modules are imported by
# bare name.
_tests = os.path.abspath(os.path.dirname(__file__))
for _p in (os.path.join(_tests, "audit"),
           os.path.join(_tests, "substrate"),
           _tests):
    sys.path.insert(0, _p)

from interference import (couplings, dark_sector, gravity, higgs_ew,
                          masses, mixing, reference, root, wiring, words)
from interference import _paths
from interference.searchlight import (embedding_uniqueness, octonions,
                                      protected_forgetting)

# Older probes use historical short spellings. These names alias the canonical
# package objects during the test session.
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
    sys.modules[_name] = _module


@pytest.fixture(scope="session")
def res():
    """Run the full derivation chain once, silently."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        c = couplings.derive(R, m)
        g = gravity.derive(R, m, mx, c)
        d = dark_sector.derive(R, g)
    return {
        "R": R,        # algebraic root (four integers + all derived rationals)
        "m": m,        # 9 masses + v_EW + bridge factor
        "mx": mx,      # CKM + PMNS
        "c": c,        # alpha_s + alpha(0) (bridge self-interference)
        "g": g,        # Newton's G + Higgs + baryogenesis + neutrinos + EW
        "d": d,        # dark sector + cosmological density bookkeeping
    }
