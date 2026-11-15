"""Mechanism probes of the couplings branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


def test_condensate_map_is_the_per_mode_reading():
    """Runs the condensate map's field map. In the real canonical norm
    the SM broken generator sum rule uses all of v_EW. A linear
    identification of the 7-singlet therefore carries v7 = v_EW
    (alpha_s at -2.5 sigma). The coded v7^2 = E[v^2] v_EW^2 is the
    per mode contraction moment of the bridge. It reproduces
    M_V = 111.648 exactly. The two contractions are 3 sigma apart.
    Every nonelectroweak scale on the ladder fails by tens of sigma.
    The look elsewhere count over bare factors is recorded: five atoms
    lie within 1 sigma, so the role fixes the choice, not the fit."""
    import contextlib
    import io

    import masses
    import root as _root
    from probes import condensate_map

    with contextlib.redirect_stdout(io.StringIO()):
        R = _root.derive()
        M = masses.derive(R)

    r = condensate_map.run(report=_silent, engine_m_t_GeV=M["m_t"] / 1e3)

    # 1. the SM sum rule: C2(2) = 3/4, three equal masses, M_W = gv/2,
    #    identical in the real and complex bookkeepings of one doublet
    sm = r["sm_real_norm"]
    assert sm["casimir_sum_over_g2v2"] == pytest.approx(0.75, abs=1e-12)
    assert sm["mw_over_gv"] == pytest.approx(0.5, abs=1e-12)
    assert (sm["real_norm_length"] / sm["complex_norm_length"]
            == pytest.approx(math.sqrt(2.0), abs=1e-12))
    # 2. the per mode contraction reproduces the coded threshold
    assert r["per_mode"]["M_V"] == pytest.approx(111.648, abs=1e-2)
    assert r["per_mode"]["M_V"] == pytest.approx(
        (math.pi / (2.0 * math.sqrt(2.0))) * r["v_EW_GeV"]
        / math.sqrt(6.0), abs=1e-9)
    assert abs(r["per_mode"]["pull"]) < 0.5
    # 3. the linear map (= the total contraction) is the excluded row
    assert r["linear_map"]["M_V"] == pytest.approx(157.894, abs=1e-2)
    assert r["linear_map"]["pull"] < -2.0
    assert r["dichotomy_separation_sigma"] == pytest.approx(2.95, abs=0.1)
    # 4. the scale menu: every nonelectroweak condensate row is excluded
    menu = r["scale_menu_pulls"]
    assert abs(menu["M_R:per_mode"]) > 50
    assert abs(menu["Lambda_conf:per_mode"]) > 50
    assert abs(r["no_threshold_pull"]) > 5
    # 5. look elsewhere count: the 1-sigma atom count is frozen
    assert len(r["atoms_within_1sigma"]) == 5
    assert "E_v2 = 1/2" in r["atoms_within_1sigma"]


@pytest.mark.slow
def test_ew_internal_reclassification():
    """G_F is internal, at ppm accuracy. The top doublet rho term is
    computed from framework outputs."""
    from probes import ew_internal
    r = ew_internal.run(report=_silent)
    assert abs(r["G_F"] / 1.1663788e-5 - 1.0) < 1e-5
    assert abs(r["d_rho_t"] - 0.00933) < 2e-4


def test_scale_ties_gluon_coupling_admits_the_radial_mode_alone(res):
    """On the engine's two loop alpha_s chain the coset threshold moves
    1/alpha_s(125 GeV) by -4.52 times the top per unit ln v. The 125 GeV
    state as the radial mode of the electroweak condensate, with G2
    broken by its own field, gives kappa_g = 1. Moving M_V gives -3.52.
    A dilaton of the transmuted scale (mu*, M_V, m_t) gives +10.31.
    Only the first lies within 2 sigma of CMS 0.92 +/- 0.08."""
    from probes import scale_ties

    r = scale_ties.run(report=_silent)
    assert r["passes"] == {"A": True, "B": False, "C": False}
    assert r["kappa"]["B"] == pytest.approx(-3.52, abs=0.02)
    assert r["kappa"]["C"] == pytest.approx(10.31, abs=0.05)
    assert r["coset_only"] / r["top"] == pytest.approx(-4.52, abs=0.01)
    assert r["alpha_s_MZ"] == pytest.approx(0.11839, abs=5e-5)
