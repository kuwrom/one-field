"""Mechanism probes of the audit branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##



## The closure dichotomy ##

def test_surplus_edges_enumerated():
    """The grammar run AWAY from the data (web completeness).  Every
    base in use is admissible.  The excluded terminus ghosts (64 and
    139 GeV) stay absent from the PDG.  The surplus channel (the free
    prediction channel) is finite and enumerated."""
    from probes import surplus_edges
    r = surplus_edges.run(report=_silent)
    assert len(r["rows"]) == 15
    assert all(g in (1, 2, 3) for g, w, m in r["falsifiers"])


def test_record_area_is_finite_and_depth_blind():
    """The Bogoliubov vacuum's entanglement per area converges with no
    cutoff (0.0227 mu). In healing units it is kappa = 0.0057 per xi^2,
    the pixel size alone (parameter_ledger.py P3 reads it at two
    chemical potentials).  Read as the closure, it places xi0 between
    1.918 and 2.93 l_Pl across the vector contact coefficients (P4)."""
    from probes import record_area

    r = record_area.run(report=_silent)
    assert r["limit"] == pytest.approx(0.0227, abs=5e-4)
    assert r["kappa"] == pytest.approx(0.00568, abs=1.5e-4)
    assert r["xi_closure"][1.0] == pytest.approx(1.918, abs=0.02)
    assert r["xi_closure"][-1.0] == pytest.approx(2.93, abs=0.03)


def test_channel_response_common_is_binary_relative_is_computed():
    """On the three component clock the common channel exchange is one
    screened term at 2 c_s, 1/(8 pi R0 xi^3) for a mu coupled quantum and
    1.418 times it for a Hartree coupled one.  The relative channels are
    two screened terms at the static stiffnesses (5 -+ sqrt 13) r/2.
    The O4 depth is 2^(3/2) times the depth in the range convention."""
    from probes import channel_response

    r = channel_response.run(report=_silent)
    assert r["eigen_error"] < 1e-8 and r["radial_residual"] < 1e-8
    assert r["hartree_over_mu"] == pytest.approx(1.41787, abs=1e-5)
    assert r["alpha_common_mu"] == pytest.approx(0.148821, abs=1e-6)
    assert r["relative_over_common"] == pytest.approx(0.12830, abs=2e-5)
    assert r["relative_ranges_xi"] == pytest.approx([1.561, 3.878], abs=2e-3)
    assert r["depth_o4"] / r["depth_range"] == pytest.approx(2**1.5, abs=1e-12)
    assert r["kappa_at_range_depth_3"] == pytest.approx(125.9, abs=0.1)
    assert r["depth_mu"] == pytest.approx(3.0, abs=1e-12)
    assert r["depth_hartree"] == pytest.approx(4.254, abs=2e-3)
