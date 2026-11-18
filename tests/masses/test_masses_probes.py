"""Mechanism probes of the masses branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The spectrum at the MTC layer, and theta's uniqueness ##

def test_mtc_spectrum():
    """The three line spectrum exists at the MTC layer, assembled
    entirely from computed data."""
    from probes import mtc_spectrum
    r = mtc_spectrum.run(report=_silent)
    assert abs(r["Q"] - 2.0 / 3.0) < 1e-12
    assert abs(r["mu_e"] - 206.77) < 0.01


def test_theta_menu_permits_only_2_9():
    """The data permit only theta = h(fund) = 2/9. Every other
    natural MTC angle leaves the positive branch."""
    from probes import theta_menu
    dev, name, mu_e, tau_e = theta_menu.run(report=_silent)
    assert dev < 1e-3
    assert abs(mu_e - 206.7683) / 206.7683 < 1e-4


## The wiring derivation ##

def test_edge_derivations():
    """The two derived mass edges. m_d is the exact on menu
    identity (aa + h10) = 83/9. The strange bridge is the norm
    conversion 32/27, unique by Schur. Its Albert gap is constructed on
    an antiassociative triple, and the isotropic share is the only one
    in the scan."""
    from probes import edge_derivations
    r = edge_derivations.run(report=_silent)
    assert r["m_d"] == "83/9" and r["bridge_sq"] == "32/27"
    assert r["delta_J_sq"] == "3/8" and r["share_hits"] == ["1/8"]
    assert r["anti_associative"] == 28


def test_koide_assignment_forced():
    """The Koide representation assignment (the insertion rule) is
    derived. Given the neutral lane premise, adjoint -> Q(c,b,t) and
    fundamental -> Q(s,c,b) are forced by two canon lemmas. The swapped
    assignment is excluded."""
    from fractions import Fraction
    from probes import koide_assignment
    r = koide_assignment.run(report=_silent)
    assert r["Q_cbt"] == Fraction(289, 432)
    assert r["Q_scb"] == Fraction(649, 972)
    assert r["neutral_stepper"] == (1, 1)
    assert r["mb_swap"] - r["mb_canon"] > 70.0     # swap ~ +13 sigma


def test_vertex_pinning_locks():
    """Each depth 4 response is pinned by its own observable. Each
    has a unique 2-sigma survivor (2, 16, 27/2) over the 29-value menu."""
    from fractions import Fraction
    from probes import vertex_pinning
    r = vertex_pinning.run(report=_silent)
    assert r["survivors"] == ([Fraction(2)], [Fraction(16)],
                              [Fraction(27, 2)])


def test_echo_naturality_regenerates():
    """The response rule. A web with zero per node depth 4
    declarations (one loop object read through channel responses)
    regenerates the canonical precision block exactly."""
    from probes import echo_naturality
    r = echo_naturality.run(report=_silent)
    assert set(r["regenerated"]) == {"inv_alpha", "m_mu", "m_tau",
                                     "v_EW"}
    assert r["responses"] == {"coupling": "2", "action": "16",
                              "amplitude": "27/2"}


def test_circulant_normalization():
    """Three claimants to the mass circulant slot, adjudicated.

    Freezes five facts.
    (a) Q = 1/3 + (B/A)^2/6 is exact and blind to theta. The Koide
        value pins the AMPLITUDE alone, and no normalisation can fix
        the phase.
    (b) The total and per mode second moments differ by the factor 2
        of the two Z3 channels (Pv2P = 2*E[v2]). This is why the web's
        E_v2 = 1/2 sits alongside Pv2P = 1.
    (c) theta = 0 forces exact two fold degeneracy at every amplitude.
        This excludes the (1,-1/2,-1/2) harmonic.
    (d) At alpha = 0 the substrate's uniform BdG circulant has
        B/A = 2(1 + 3r/2)/(1 - r) with r = g1/g0 (14.07 at
        r = 1/sqrt 2, and >= 2 for every r in [0, 1)). There the
        positive gap branch and the modulational instability are
        DISJOINT.
    (e) A complex Z3 phase is admissible. It reaches (sqrt 2, 2/9) at a
        SINGLE physical point, which is modulationally stable (full
        6x6 BdG, PART F; direct dynamics arbitration, PART G). Its 'two
        inverse images' are one theory under the exact redundancy
        (g1, alpha) -> (-g1, alpha - pi). The disjointness holds for
        every admissible (r, alpha), and the canonical point is a
        stable uniform phase carrying the Koide circulant.
    A per channel pairing reading reports coexistence instead. That is
    an artifact of the pairing, valid only at alpha = 0. The effective
    coupling route has no implementation."""
    from probes import circulant_normalization as cn
    r = cn.run(report=_silent)
    assert r["identity_max_dev"] < 1e-12
    assert abs(r["Pv2P_total"] - 1.0) < 1e-12
    assert abs(r["E_v2_mode"] - 0.5) < 1e-12
    assert r["theta0_degenerate"] is True
    assert abs(r["Q_canon"] - 2.0 / 3.0) < 1e-12
    assert abs(r["Q_harmonic"] - 0.5) < 1e-12
    assert r["BdG_BA"] >= 2.0
    assert r["r_needed_for_sqrt2"] < 0.0
    assert r["gap_min"] > 0.0
    assert r["reconcilable"] is False
    # PART F + G: the canonical point is reached and is stable. Its
    # 'two inverse images' are ONE theory (exact redundancy, asserted).
    # The dynamics arbitration confirms this: the canonical seeded mode
    # stays flat, and the alpha = 0 control grows at g1 rho0. Frozen as
    # a consequence of the hermitian equation.
    assert abs(r["BA_alpha"] - math.sqrt(2.0)) < 1e-9
    assert abs(r["theta_alpha"] - 2.0 / 9.0) < 1e-9
    assert r["alpha_point_unstable"] is False
    assert r["alpha_images_all_stable"] is True
    assert r["alpha_point_positive_gaps"] is True
    assert r["inverse_image_redundancy_dev"] < 1e-9
    assert abs(r["dynamics_rate_canonical"]) < 0.02
    assert 0.5 < r["dynamics_rate_control"] < 0.85
    # the closed form circulant must track nls_soliton's live matrix
    assert cn.test_matches_live_substrate() is True


def test_harmonic_single_source():
    """The canonical harmonic (B, theta) is defined in ONE place:
    root.py. protected_forgetting.py carries its own copies of the
    values so it runs without the engine. This test checks that the
    copies cannot drift. They equal root's canonical values, and the
    module carries a drift assertion at module level, which runs when
    the module is imported with the engine on the path."""
    import root
    import protected_forgetting as pf
    assert pf.B_CANON == root.BA_ratio
    assert pf.THETA_CANON == float(root.h10)
    v = pf.harmonic()
    # the electron gap that anchors M_Pl is the k=1 component
    assert abs((1.0 + v[1]) - root._Delta_e) < 1e-15
