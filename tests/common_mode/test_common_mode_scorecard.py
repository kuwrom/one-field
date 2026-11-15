"""Scorecard tests of the common mode branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Quark masses  (6 predictions)
# ═══════════════════════════════════════════════════════════════════════



# ═══════════════════════════════════════════════════════════════════════
#  Higgs running coordinate
# ═══════════════════════════════════════════════════════════════════════

def test_higgs_transport_coordinate(res):
    """The candidate boundary is transported in one MS-bar scheme."""
    g = res["g"]
    assert g["lambda_mt_transport"] == pytest.approx(0.1265078, abs=2e-7)
    assert g["higgs_transport_coordinate_GeV"] == pytest.approx(
        123.8501, abs=2e-3)


# ═══════════════════════════════════════════════════════════════════════
#  Gravity sector  (Newton constant from the (7,26) bridge heat kernel)
# ═══════════════════════════════════════════════════════════════════════

def test_newton_internal_reading(res):
    """Freeze the gravity normalisation's internal reading under the
    diagnostic cutoff t0 = d11 l_Pl^2 (FACT 7; standard per state
    weights). Not an independent G test."""
    assert abs(res["g"]["G_ratio_UV"] - 0.698131701) < 1e-8


def test_newton_phase_independence(res):
    """Proca = massless vector + eaten scalar: the symmetric and broken
    ledgers are identical."""
    assert abs(res["g"]["G_ratio_broken"] - res["g"]["G_ratio_UV"]) < 1e-12


# ═══════════════════════════════════════════════════════════════════════
#  Neutrino sector (structural rank bound; ordering not selected)
# ═══════════════════════════════════════════════════════════════════════

def test_neutrino_rank_bound(res):
    """Minimal type-I with two identified RH neutrinos has rank <= 2."""
    g = res["g"]
    assert g["neutrino_mass_mechanism_assumption"] == \
        "minimal type-I seesaw only"
    assert g["n_RH"] == 2
    assert g["light_mass_rank_max"] == 2
    assert g["massless_light_states_min"] == 1
    assert g["massless_state_by_ordering"] == {
        "normal": "m1", "inverted": "m3"}


def test_neutrino_ordering_is_not_selected(res):
    """The rank bound permits both normal and inverted ordering."""
    assert res["g"]["ordering"] == "undetermined"


# ═══════════════════════════════════════════════════════════════════════
#  Dark sector  (Omega_DM / Omega_b = 2 pi - 1)
# ═══════════════════════════════════════════════════════════════════════

def test_dm_baryon_ratio(res):
    """Omega_DM / Omega_b = 2 pi - 1 (bridge self-interference)."""
    from interference import wiring

    DM_b = res["d"]["DM_baryon_ratio"]
    expected = 2.0 * math.pi - 1.0
    assert math.isclose(DM_b, expected, rel_tol=1e-12)
    assert res["d"]["f_baryonic"] == wiring.normalized_current_response(1.0)
    # Within 2 sigma of Planck
    assert abs(res["d"]["pull_DM_baryon"]) < 2.0


def test_xi_bridge_is_one_over_48_pi(res):
    """The label's coupling xi_bridge = alpha_G_2 * E[v^2] * h_bridge = 1/(48 pi).
    Only the labelled bridge reads it.  The amount does not."""
    from interference import boundary
    assert boundary.xi_of("amount", 1.0) == 0.0 and boundary.xi_of("boundary", 1.0) == 1.0
    R = res["R"]
    assert math.isclose(R["xi_bridge"], 1.0 / (48.0 * math.pi), abs_tol=1e-14)
    assert math.isclose(R["E_v2"], 0.5, abs_tol=1e-12)


def test_newton_constant_within_3sigma():
    """G = (3Q_l-1) hbar c/M_b^2, the anchor inversion read through the
    common mode (identification): within 3 sigma of CODATA 2022
    (-2.43 sigma)."""
    from root import G_PRED, G_CODATA, G_REL_SIGMA
    assert abs(G_PRED / G_CODATA - 1.0) < 3 * G_REL_SIGMA


def test_neutrino_spectrum(res):
    """Under the normal ordering benchmark, m1 = 0 plus measured
    splittings gives the displayed sum, m_beta, and m_bb band."""
    import root as _r
    pmns = res["mx"]
    s12sq = pmns.get('sin2_12', pmns.get('pmns_sin2_12'))
    s13sq = pmns.get('sin2_13', pmns.get('pmns_sin2_13'))
    m2 = math.sqrt(_r.NUFIT_OSCILLATION['dm2_21'][0])
    m3 = math.sqrt(_r.NUFIT_OSCILLATION['dm2_31'][0])
    c13sq = 1.0 - s13sq
    m_beta = math.sqrt(c13sq * s12sq * m2**2 + s13sq * m3**2)
    t2, t3 = s12sq * c13sq * m2, s13sq * m3
    assert abs((m2 + m3) * 1e3 - 58.78) < 0.05
    assert abs(m_beta * 1e3 - 8.81) < 0.05
    assert abs(abs(t2 - t3) * 1e3 - 1.51) < 0.05
    assert abs((t2 + t3) * 1e3 - 3.70) < 0.05


def test_alpha_pins_higgs_coordinate_not_sigma(res):
    """The alpha_G2(M_Pl) pins stated in root.py, measured.
    (1) the transported Higgs coordinate responds to the bridge boundary.
    (2) the amount is minimally coupled (xi = 0), so Sigma = 27 does not
    depend on alpha."""
    import numpy as np
    import higgs_ew as hew
    import root as rt

    m, c, g = res["m"], res["c"], res["g"]
    v_EW = m["v_EW_pred_GeV"]
    mt = m["m_t"] / 1000.0
    log_mt_mz = math.log(mt / g["M_Z_derived"])
    g1_mz = math.sqrt(4 * math.pi * g["alpha_em_MZ"]
                      / (1 - g["sin2W"])) * math.sqrt(5.0 / 3.0)
    g2_mz = math.sqrt(4 * math.pi * g["alpha_em_MZ"] / g["sin2W"])
    g1 = hew._run_gauge_one_loop(g1_mz, 41.0/10.0, log_mt_mz)
    g2 = hew._run_gauge_one_loop(g2_mz, -19.0/6.0, log_mt_mz)
    g3 = math.sqrt(4 * math.pi * c["alpha_s_mt"])
    yt = g["yt_mt_transport"]
    t_pl = math.log(rt.M_Pl_GeV / mt)

    def m_higgs(lam_pl, n=8000):
        lam_it = 0.0
        for _ in range(3):
            up = hew._run_rge(np.array([g1, g2, g3, yt, lam_it]),
                              0.0, t_pl, hew._beta_sm, n)
            dn = hew._run_rge(np.array([up[0], up[1], up[2], up[3],
                                        lam_pl]),
                              t_pl, 0.0, hew._beta_sm, n)
            lam_it = dn[4]
        return math.sqrt(2.0 * abs(lam_it)) * v_EW

    m0 = g["higgs_transport_coordinate_GeV"]
    d_plus = m_higgs(g["lam_Pl_boundary"] * 1.21) - m0
    d_minus = m_higgs(g["lam_Pl_boundary"] * 0.81) - m0
    assert d_plus < -0.5
    assert d_minus > 0.5

    # (2) the stiffness depends on no coupling: the amount is minimally
    # coupled (boundary.xi_of), so Sigma is the count 162/6 = 27 at any alpha
    assert abs(g["Sigma_total_UV"] - 27.0) < 1e-12
    assert rt.WEB["a1_bridge"].terms == []
    # executable: the counted list at alpha and at 1.1 alpha gives one Sigma
    from interference import boundary
    booked = boundary.rows(rt.d10, rt.d11, rt.n7, rt.n26, scalars="light")
    xi = float(rt.xi_bridge)
    s0 = boundary.sigma(booked, 0.0, 1.0, xi)
    s1 = boundary.sigma(booked, 0.0, 1.0, 1.1 * xi)
    assert s0 == s1 == pytest.approx(27.0, abs=1e-12)


def test_two_F4_singlets_available_for_RH_identification():
    """The 26 decomposition contains two singlets. Identifying them as
    RH neutrinos gives a rank at most two matrix only in minimal type-I."""
    import root

    assert root.n26 - 3 * (root.d11**2 - 1) == 2


def test_no_higgs_mass_leaves_the_chain(res):
    """The mixed scheme 125.2965 reading is not canonical: no Higgs
    mass prediction leaves the derivation."""
    assert "mH_pred" not in res["g"]
