"""The probes as tests: every mechanism claim is an assertion.

The scorecard suite (test_interference.py) freezes the NUMBERS. This
suite freezes the MECHANISMS.  Each probe runs on the theory's own
equations and its registered conclusion is asserted here, so a change
that silently breaks a mechanism claim fails CI exactly like a change
that moves a mass.

The substrate-dynamics probes evolve the NLS equation for
minutes and carry the `slow` marker. The default run skips them:

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
    """E8 > G2(1) x F4(1) at k = 3 is the unique zero-residual wiring.
    The dictionary is read off the winner, never input."""
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
    chirality blocks. The windingless mass is gapped.  Handedness is
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


# ── Statistics: the level counted, the winding carried ──────────────

def test_bookkeeping_level_and_theta():
    """Child anomaly vanishes. WZ level = 3 by counting. The lightest
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

def test_effective_action_stabilizer():
    """DP #4 branch (i), closed constructively: the exact CP^1 split
    and Mermin-Ho relation hold on the repo's texture (verified with
    resolution convergence), and the induced Faddeev-Skyrme quartic
    (coefficient rho/8m^2, m the measured ring-down gap) produces a
    Derrick minimum at the healing scale, with the local limit
    self-consistent at its own minimum."""
    from probes import effective_action
    r = effective_action.run(report=_silent)
    assert r["R_star_m"] > 1.0
    assert 0.5 < r["R_star_xi"] < 2.0


def test_surplus_edges_enumerated():
    """The grammar run AWAY from the data (web completeness): every
    used base is admissible, the excluded terminus ghosts (64, 139
    GeV) stay absent from the PDG, and the surplus channel is finite
    and enumerated (the free prediction channel)."""
    from probes import surplus_edges
    r = surplus_edges.run(report=_silent)
    assert len(r["rows"]) == 15
    assert all(g in (1, 2, 3) for g, w, m in r["falsifiers"])


def test_selector_theorem_binary():
    """The idempotent-measure space on Z3 is exactly {Haar, delta_e}:
    the confinement dichotomy as a two-element theorem."""
    from probes import selector_theorem
    found = selector_theorem.run(report=_silent)
    assert len(found) == 2


# ── Substrate dynamics (the classical-1D arc. Minutes of PDE) ───────

@pytest.mark.slow
def test_knot_formation_and_control():
    """A winding seed grows and localizes. With g1 = 0 the same seed
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
    """Standing character is universal. A unique attractor is NOT
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
    its docstring. The fringe question is interpretive)."""
    from probes import spectrum
    spectrum.run(report=_silent)


@pytest.mark.slow
def test_ew_internal_reclassification():
    """G_F internal at ppm. The top-doublet rho term computed from
    framework outputs."""
    from probes import ew_internal
    r = ew_internal.run(report=_silent)
    assert abs(r["G_F"] / 1.1663788e-5 - 1.0) < 1e-5
    assert abs(r["d_rho_t"] - 0.00933) < 2e-4


@pytest.mark.slow
def test_knot_charge_not_protected():
    """pi_3 winding is NOT a substrate-level invariant: the estimator
    converges, the saturated knot unwinds through amplitude-node
    events, and no Noether charge backs internal rotation (the Q-lump
    loophole is closed).  Frozen witness for DERIVATION_PROGRAMS #4."""
    from probes import knot_charge
    knot_charge.main()  # four assertions live inside


@pytest.mark.slow
def test_skyrmion_3d_formation():
    """3D knot formation (the review's sim3d experiment as a house
    probe): the seeded texture grows only under the true G2 coupling
    (control disperses), the growth rate matches the computed BdG
    instability, and the common mode dips (the shadow).  Charge
    protection is NOT claimed. knot_charge.py measures the
    unwinding. Protection belongs to the confined phase."""
    from probes import skyrmion_3d
    r = skyrmion_3d.run(report=_silent)   # frozen assertions inside
    assert r["growth"] > 5.0 and r["decay"] < 0.05
    assert (0.5 * r["sigma_bdg"] < r["sigma_measured"]
            <= 1.1 * r["sigma_bdg"])


@pytest.mark.slow
def test_confined_phase_entry():
    """The layer-bridge entry test (SUBSTRATE_CONJECTURES #1): the Z3
    disorder operator measured exactly (plaquette vortex density).
    The imprinted ring registers as the positive control. The
    classical endpoint is defect-free (perimeter law), so the
    confined phase is not a classical-field property.  The
    conjecture itself stays OPEN."""
    from probes import confined_phase
    confined_phase.run(report=_silent)


@pytest.mark.slow
def test_knot_spectrum_tower():
    """The saturated lump's ring-down tower is gapped, with the first
    internal excitation at the relative-sector binding scale: no
    surplus states below the lock-in scale (feeds the D4 prediction
    channel)."""
    from probes import knot_spectrum
    r = knot_spectrum.run(report=_silent)
    assert len(r["pk_rel"]) >= 1
