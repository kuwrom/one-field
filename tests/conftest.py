"""Shared fixtures for the scorecard test suite.

The full derivation chain is run once per test session, silently.  Every
test then asserts against the cached results, so the suite is fast and
the layer-by-layer console output never pollutes pytest.

Run with:
    pytest -q
"""

import contextlib
import io
import os
import sys

import pytest

# Add the sibling interference/ package to sys.path so its modules
# (root, masses, ...) import as bare names, matching how run.py
# bootstraps the chain.
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))
import _paths  # noqa: E402,F401  (adds interference/ and searchlight/)
# Add tests/ itself so the mechanism probes import as `probes.*`, plus
# the two scrutiny directories: substrate/ (the NLS conjecture) and
# audit/ (grammar + closure search), whose modules are imported by
# bare name.
_tests = os.path.abspath(os.path.dirname(__file__))
for _p in (os.path.join(_tests, "audit"),
           os.path.join(_tests, "substrate"),
           _tests):
    sys.path.insert(0, _p)

import root
import masses
import mixing
import couplings
import gravity
import dark_sector


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
        "d": d,        # dark sector + CC scale
    }
