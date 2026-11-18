"""Scorecard tests of the couplings branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Strong coupling  (1 prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_alpha_s_at_MZ(res):
    """alpha_s(M_Z) from the embedding index chain SU(3) C G_2 C E_8."""
    alpha_s_MZ = res["c"]["alpha_s_MZ_thresh"]
    pdg = 0.1180
    assert err_pct(alpha_s_MZ, pdg) < 0.76      # within 1 sigma, 0.0009


def test_full_record_strong_cp_closure(res):
    """The production chain exposes the complete vacuum character,
    not separately basis dependent bare theta and mass phases."""
    sc = res["c"]["strong_cp"]
    assert res["c"]["bar_theta"] == 0.0
    assert sc["physical_angle"] == "bar_theta"
    assert sc["phase_rule_status"] == "CONSTITUTIVE ONTOLOGY"
    assert sc["vacuum_phase_automorphisms"] == (1.0 + 0j,)
    assert all(value == 1.0 + 0j
               for value in sc["winding_character_samples"].values())
    assert not sc["theta_QCD_separately_defined"]
    assert not sc["mass_determinant_phase_separately_defined"]


def test_smdr_beta_reference_point():
    """Freeze the yb=ytau=0 port against SMDR v1.3 betas.c."""
    import higgs_ew

    y = np.array([0.46232945922286817, 0.6470936480109292,
                  1.1649178681893946, 0.9350096081781627,
                  0.1266556161139303])
    got = higgs_ew._beta_sm(y)
    expected = np.array([
        0.0020276825305886822 / math.sqrt(3.0/5.0),
        -0.0052288603702881362,
        -0.0722744362694192483,
        -0.0537308630890171654,
        -0.0199725406547064908,
    ])
    assert np.allclose(got, expected, rtol=0.0, atol=5e-16)


# ═══════════════════════════════════════════════════════════════════════
#  Electroweak (derived)  (sin^2 theta_W, M_Z, M_W)
# ═══════════════════════════════════════════════════════════════════════

def test_sin2_theta_W(res):
    """sin^2 theta_W = 3/13 + Q₀²·α/(2π), depth 1 h₇ + depth 2 h₇/d₁₁²."""
    sin2W = res["g"]["sin2W"]
    # canonical: base 3/13 + depth 1 h₇·α/(2π) + depth 2 (h₇/d₁₁²)·α/(2π)
    # total echo coefficient h₇(1+1/d₁₁²) = Q₀² = 4/9
    import root as _r
    expected = 3.0/13.0 + (4.0/9.0) * _r.alpha_phys / (2.0 * math.pi)
    assert math.isclose(sin2W, expected, rel_tol=1e-12)
    assert abs(sin2W - 0.23129) / 4e-5 < 2.0      # within 2 sigma of PDG 2024 global fit
    # PDG: 0.23122
    assert err_pct(sin2W, 0.23129) < 0.5


def test_M_Z_derived(res):
    """M_Z from the declared MS-bar transport imports: within 0.02%."""
    M_Z = res["g"]["M_Z_derived"]
    M_Z_PDG = res["g"]["M_Z_PDG"]
    assert err_pct(M_Z, M_Z_PDG) < 0.02


def test_M_W_derived(res):
    """M_W from alpha(0), v_EW, and the declared dr_hat_W import:
    within 2 sigma of the 80.3625(77) world average (PDG 2026)."""
    M_W = res["g"]["M_W_derived"]
    M_W_PDG = res["g"]["M_W_PDG"]
    assert abs(M_W - M_W_PDG) / res["g"]["M_W_err"] < 2.0
    assert (M_W - M_W_PDG) / res["g"]["M_W_err"] == pytest.approx(-0.86, abs=0.01)
    assert err_pct(M_W, M_W_PDG) < 0.5


def test_singh_ratio_is_16(res):
    """alpha_s / alpha_em = (8/3) * C_2(26) = (8/3) * 6 = 16."""
    import root
    d10, d11 = root.d10, root.d11
    assert root.charge_trace == Fraction(d10**2 + 1 + d11, d11) == Fraction(8, 3)
    assert root.singh_ratio == 16
    assert root.alpha_EM == pytest.approx(root.alpha_G2_WZW / 16, rel=1e-15)


# ═══════════════════════════════════════════════════════════════════════
#  alpha(0) bridge self-interference  (h=1, D^2=1, c=0)
# ═══════════════════════════════════════════════════════════════════════

def test_alpha_bridge_chain_is_one_loop_exact(res):
    """alpha(0) bridge derivation: h_bridge = 1, D^2 = 1, c_coset = 0."""
    bridge = res["c"]["bridge"]
    assert math.isclose(bridge["h_bridge"], 1.0, abs_tol=1e-12)
    assert math.isclose(bridge["h_G2_7"], 2 / 5, abs_tol=1e-12)
    assert math.isclose(bridge["h_F4_26"], 3 / 5, abs_tol=1e-12)
    assert math.isclose(bridge["D2_local"], 1.0, abs_tol=1e-12)
    assert math.isclose(bridge["g_bridge_sq"], 1.0, abs_tol=1e-12)
    assert math.isclose(bridge["c_coset"], 0.0, abs_tol=1e-12)
    # Depth-1 closed form 256(2 pi - 1)/pi^2 is the truncation. The full
    # value includes the FORCED depth 4 reading of the e<->q carrier and
    # solves the cubic  x^3 = (512/pi)[(1 - 1/(2 pi)) x^2 - 1/(2 pi^2)].
    x = bridge["inv_alpha_phys"]
    depth1 = 256.0 * (2.0 * math.pi - 1.0) / math.pi ** 2
    assert x < depth1  # the depth 4 response lowers 1/alpha
    cubic = x**3 - (2**9 / math.pi) * ((1.0 - 1.0/(2.0*math.pi)) * x**2
                                       - 1.0/(2.0*math.pi**2))
    assert abs(cubic) < 1e-6
    # matches the Berkeley Cs recoil measurement 137.035999046(27)
    assert abs(x - 137.035999046) < 27e-9


def test_vew_matches_GF_at_1e7(res):
    """v_EW (with the FORCED 16-term and the m_e anchor) matches the
    G_F determination at the 1e-7 level."""
    v_target = (math.sqrt(2.0) * 1.1663787e-5) ** -0.5
    assert abs(res["m"]["v_EW_pred_GeV"] / v_target - 1.0) < 5e-7


def test_scale_landmarks(res):
    """The scale ladder, with M_V the connector's per mode reading."""
    import root as _r
    M_V = res["c"]["M_V_GeV"]
    assert abs(M_V - 111.648) < 1e-2
    g_G2 = math.sqrt(4.0 * math.pi * res["c"]["alpha_s_mu_star"])
    n_broken = 2 * _r.n7 - (_r.d11**2 - 1)
    v_break = res["m"]["v_EW_pred_GeV"] / math.sqrt(_r.d10)
    assert M_V == pytest.approx(
        g_G2 * v_break * math.sqrt(_r.C2_7 / n_broken), abs=1e-12)
    M_G2 = _r.M_Pl_GeV * math.exp(-2.0 * math.pi**2)
    assert abs(M_G2 / 3.27e10 - 1.0) < 0.01


def test_condensate_map_is_the_per_mode_reading(res):
    """The G2 breaking scale is the connector's PER MODE reading of the
    electroweak condensate: v7^2 = E[v^2] v_EW^2 with E[v^2] = 1/2, so
    M_V = g v_EW sqrt(C2(7)/6) with g = pi/(2 sqrt(2))."""
    import root

    v_ew = root.M_Pl_GeV * math.exp(-root.WEB.state["S_quark"])
    m_v = (math.pi / (2.0 * math.sqrt(2.0))) * v_ew * math.sqrt(
        float(root.E_v2) * float(root.C2_7) / 6.0)
    assert res["c"]["M_V_GeV"] == pytest.approx(m_v, abs=1e-2)
    assert float(root.E_v2) == 0.5
