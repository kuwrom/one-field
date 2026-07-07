"""The probes as tests: every mechanism claim is an assertion.

The scorecard suite (test_interference.py) freezes the NUMBERS; this
suite freezes the MECHANISMS.  Each probe runs on the theory's own
equations and its registered conclusion is asserted here, so a change
that silently breaks a mechanism claim fails CI exactly like a change
that moves a mass.

The five substrate-dynamics probes evolve the NLS equation for
minutes and carry the `slow` marker; the default run skips them:

    pytest -q              # fast probes included
    pytest -q -m slow      # the substrate-dynamics arc
"""

import math
import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))


def _silent(*args, **kwargs):
    pass


# ── The wiring and the dictionary ────────────────────────────────────

def test_wiring_found_not_stated():
    """E8 > G2(1) x F4(1) at k = 3 is the unique zero-residual wiring;
    the dictionary is read off the winner, never input."""
    from probes import wiring_scan
    r = wiring_scan.run(report=_silent)
    assert r["lattice"] == "E8"
    assert r["pair"] == "G2(1) x F4(1)"
    assert r["k"] == 3
    assert (r["n7"], r["n26"]) == (7, 26)
    assert abs(r["d10"] - 2.0) < 1e-9 and abs(r["d11"] - 3.0) < 1e-9


# ── Chirality: V-A computed on the knot ──────────────────────────────

def test_zero_mode_chirality():
    """Winding +1 and -1 bind one zero mode each, in opposite
    chirality blocks; the windingless mass is gapped.  Handedness is
    the winding sign, the propagated Fano bit."""
    from probes import zero_mode
    r = zero_mode.run(report=_silent)
    plus = r["winding +1 (chiral mass)"]
    minus = r["winding -1 (reversed)"]
    zero = r["winding 0 (sigma_z mass)"]
    assert plus["zero_mode"] and minus["zero_mode"]
    assert plus["gamma"] * minus["gamma"] < 0
    assert not zero["zero_mode"]


def test_orientation_bit():
    """One orientation datum: the CP-odd sign flips coherently with
    the Fano orientation."""
    from probes import orientation_bit
    r = orientation_bit.run(report=_silent)
    assert r["flip_coherent"]


# ── Statistics: the level counted, the texture protected ────────────

def test_bookkeeping_level_and_theta():
    """Child anomaly vanishes; WZ level = 3 by counting; the lightest
    Z3 sector-changer is h = 2/9 (internal assertions)."""
    from probes import bookkeeping
    bookkeeping.run(report=_silent)


def test_winding_texture_unit_charge():
    """The 3D texture carries unit topological charge (pi_3)."""
    from probes import winding_texture
    B = winding_texture.run(report=_silent)
    assert abs(abs(B) - 1.0) < 0.15


# ── The spectrum at the MTC layer, and theta's uniqueness ───────────

def test_mtc_spectrum():
    """The three-line spectrum exists at the MTC layer, assembled
    entirely from computed data."""
    from probes import mtc_spectrum
    r = mtc_spectrum.run(report=_silent)
    assert abs(r["Q"] - 2.0 / 3.0) < 1e-12
    assert abs(r["mu_e"] - 206.77) < 0.01


def test_theta_menu_permits_only_2_9():
    """The data permits only theta = h(fund) = 2/9: every other
    natural MTC angle leaves the positive branch."""
    from probes import theta_menu
    dev, name, mu_e, tau_e = theta_menu.run(report=_silent)
    assert dev < 1e-3
    assert abs(mu_e - 206.7683) / 206.7683 < 1e-4


# ── The selector dichotomy ───────────────────────────────────────────

def test_selector_theorem_binary():
    """The idempotent-measure space on Z3 is exactly {Haar, delta_e}:
    the confinement dichotomy as a two-element theorem."""
    from probes import selector_theorem
    found = selector_theorem.run(report=_silent)
    assert len(found) == 2


# ── Substrate dynamics (the classical-1D arc; minutes of PDE) ───────

@pytest.mark.slow
def test_knot_formation_and_control():
    """A winding seed grows and localizes; with g1 = 0 the same seed
    disperses.  Matter forms because of the G2 coupling."""
    from probes import knot
    r = knot.run(report=_silent)
    assert r["growth_true"] > 5.0
    assert r["growth_ctrl"] < 1.0


@pytest.mark.slow
def test_stationary_persistence():
    """The saturated lump is alive at t = 40 with a sub-scale internal
    frequency (the first `mass is rotation rate` measurement)."""
    from probes import stationary
    r = stationary.run(report=_silent)
    assert r["persistent"]
    assert 0.0 < r["omega_rel"] < r["gap_scale"]


@pytest.mark.slow
def test_universality_registered_outcome():
    """Standing character is universal; a unique attractor is NOT
    established in 1D (the registered outcome)."""
    from probes import universality
    same = universality.run(report=_silent)
    assert same is False


@pytest.mark.slow
def test_dynamics_erosion_and_shadow():
    """Bias erosion and shadow formation, from the substrate equation."""
    from probes import dynamics
    r = dynamics.run(report=_silent)
    assert r["erosion"] and r["shadow"] and r["localized"]


@pytest.mark.slow
def test_spectrum_runs():
    """The spectroscopy probe evolves stably (meanings pre-stated in
    its docstring; the fringe question is interpretive)."""
    from probes import spectrum
    spectrum.run(report=_silent)


@pytest.mark.slow
def test_ew_internal_reclassification():
    """G_F internal at ppm; the top-doublet rho term computed from
    framework outputs."""
    from probes import ew_internal
    r = ew_internal.run(report=_silent)
    assert abs(r["G_F"] / 1.1663788e-5 - 1.0) < 1e-5
    assert abs(r["d_rho_t"] - 0.00933) < 2e-4
