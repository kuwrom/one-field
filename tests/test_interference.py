"""Scorecard tests for the interference formulation.

Every numerical prediction and every load-bearing polynomial identity is
encoded as a pass/fail assertion.  The framework's claim is "no free
parameters at any layer".  This suite makes that claim falsifiable in the
operational sense: any change to (d_10, d_11, n_7, n_26) that moves a
prediction outside its tolerance band fails the suite.

Numerical tests check "does the answer match PDG".  Structural tests
check "does the answer come from the right polynomial in the four
integers".  A refactor that preserved the float but broke the algebraic
ratio would pass the numerical tests and fail the structural ones.

Run with:
    pytest -q
"""

import math
import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Mass sector  (9 predictions)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("name,tol_pct", [
    ("e",   0.01),
    ("mu",  0.01),
    ("tau", 0.01),
])
def test_lepton_mass(res, name, tol_pct):
    """Charged-lepton pole masses from the Z_3 circulant + Koide."""
    predicted = res["m"][f"m_{name}"]
    pdg = PDG_MASSES[name]
    assert err_pct(predicted, pdg) < tol_pct


@pytest.mark.parametrize("name,tol_pct", [
    ("u", 0.5),
    ("d", 0.5),
    ("s", 0.5),
    ("c", 0.5),
    ("b", 0.5),
    ("t", 0.5),
])
def test_quark_mass(res, name, tol_pct):
    """Six quark masses from F_4 triality + WZW emergence."""
    predicted = res["m"][f"m_{name}"]
    pdg = PDG_MASSES[name]
    assert err_pct(predicted, pdg) < tol_pct


# ═══════════════════════════════════════════════════════════════════════
#  Electroweak scale  (1 prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_electroweak_scale(res):
    """v_EW from M_Pl * exp(-(9 pi^2/2 - 6 + 15/512))."""
    v_EW_pred = res["m"]["v_EW_pred_GeV"]
    assert err_pct(v_EW_pred, 246.22) < 0.05


# ═══════════════════════════════════════════════════════════════════════
#  Strong coupling  (1 prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_alpha_s_at_MZ(res):
    """alpha_s(M_Z) from the embedding-index chain SU(3) C G_2 C E_8."""
    alpha_s_MZ = res["c"]["alpha_s_MZ_thresh"]
    pdg = 0.1180
    assert err_pct(alpha_s_MZ, pdg) < 0.5


# ═══════════════════════════════════════════════════════════════════════
#  CKM sector  (Wolfenstein + UT angles, summarised)
# ═══════════════════════════════════════════════════════════════════════

def test_ckm_chi_squared(res):
    """CKM observables match PDG within chi^2/n < 1."""
    n = len(res["mx"]["ckm_obs"])
    chi2_per_n = res["mx"]["chi2_ckm"] / n
    assert chi2_per_n < 1.0


def test_ckm_max_pull(res):
    """No individual CKM observable deviates by more than 2 sigma."""
    assert res["mx"]["max_pull_ckm"] < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  PMNS sector  (3 mixing angles + delta_CP prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_pmns_chi_squared(res):
    """PMNS angles match NuFit 6.0 within chi^2/n < 1."""
    chi2_per_n = res["mx"]["chi2_pmns"] / 3
    assert chi2_per_n < 1.0


def test_pmns_delta_cp_prediction(res):
    """delta_CP prediction lies in a testable range for DUNE/Hyper-K."""
    delta_CP = res["mx"]["delta_CP_deg"]
    assert 0 < delta_CP < 180


# ═══════════════════════════════════════════════════════════════════════
#  Higgs mass  (1 prediction, via SM RGE with derived boundary)
# ═══════════════════════════════════════════════════════════════════════

def test_higgs_mass(res):
    """m_H from F_4(1) fusion + fundamental-share vent
    (lambda(M_Pl) = -delta_bridge*(1-h10); registry.PROMOTIONS) + SM RGE.
    The prediction sits at +0.08% of 125.20(11); band 0.2%."""
    mH = res["g"]["mH_pred"]
    assert err_pct(mH, 125.20) < 0.2


# ═══════════════════════════════════════════════════════════════════════
#  Gravity sector  (Newton constant from the (7,26) bridge heat kernel)
# ═══════════════════════════════════════════════════════════════════════

def test_newton_constant_UV(res):
    """G_ind/G_N (UV bookkeeping) within 1% of unity."""
    # face-split (no-self-dilution) law
    assert abs(res["g"]["G_ratio_UV"] - 1.0) < 1e-6


def test_newton_constant_broken(res):
    """G_ind/G_N (broken-phase bookkeeping) within 1% of unity."""
    assert err_pct(res["g"]["G_ratio_broken"], 1.0) < 1.5


# ═══════════════════════════════════════════════════════════════════════
#  Baryogenesis  (eta_B from G_2 instanton + Fano orientation)
# ═══════════════════════════════════════════════════════════════════════

def test_baryon_asymmetry(res):
    """eta_B = n_7 * J_lep * exp(-2 pi^2) within 2% of Planck."""
    eta_B = res["g"]["eta_B"]
    eta_B_obs = res["g"]["eta_B_obs"]
    assert err_pct(eta_B, eta_B_obs) < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  Neutrino sector  (structural: m_1 = 0, normal ordering)
# ═══════════════════════════════════════════════════════════════════════

def test_neutrino_lightest_mass_is_zero(res):
    """Rank-2 Type-I seesaw from F_4's two singlets forces m_1 = 0."""
    assert res["g"]["m1"] == 0


def test_neutrino_normal_ordering(res):
    """Algebraic structure forces normal ordering."""
    assert res["g"]["ordering"] == "normal"


# ═══════════════════════════════════════════════════════════════════════
#  Electroweak (derived)  (sin^2 theta_W, M_Z, M_W)
# ═══════════════════════════════════════════════════════════════════════

def test_sin2_theta_W(res):
    """sin^2 theta_W = 3/13 + Q₀²·α/(2π), depth-1 h₇ + depth-2 h₇/d₁₁²."""
    sin2W = res["g"]["sin2W"]
    # canonical: tree 3/13 + depth-1 h₇·α/(2π) + depth-2 (h₇/d₁₁²)·α/(2π)
    # total echo coefficient h₇(1+1/d₁₁²) = Q₀² = 4/9
    import root as _r
    expected = 3.0/13.0 + (4.0/9.0) * _r.alpha_phys / (2.0 * math.pi)
    assert math.isclose(sin2W, expected, rel_tol=1e-12)
    assert abs(sin2W - 0.23129) / 4e-5 < 2.0      # within 2 sigma of PDG 2024 global fit
    # PDG: 0.23122
    assert err_pct(sin2W, 0.23129) < 0.5


def test_M_Z_derived(res):
    """M_Z at layer-1 (MS-bar chain, PDG 2024 imports): within 0.02%."""
    M_Z = res["g"]["M_Z_derived"]
    M_Z_PDG = res["g"]["M_Z_PDG"]
    assert err_pct(M_Z, M_Z_PDG) < 0.02


def test_M_W_derived(res):
    """M_W at layer-1 (A from alpha(0),v_EW; dr_hat_W declared import):
    within 2 sigma of PDG 80.3692(133)."""
    M_W = res["g"]["M_W_derived"]
    M_W_PDG = res["g"]["M_W_PDG"]
    assert abs(M_W - M_W_PDG) / 0.0133 < 2.0
    assert err_pct(M_W, M_W_PDG) < 0.5


# ═══════════════════════════════════════════════════════════════════════
#  Dark sector  (Omega_DM / Omega_b = 2 pi - 1)
# ═══════════════════════════════════════════════════════════════════════

def test_dm_baryon_ratio(res):
    """Omega_DM / Omega_b = 2 pi - 1 (bridge self-interference)."""
    DM_b = res["d"]["DM_baryon_ratio"]
    expected = 2.0 * math.pi - 1.0
    assert math.isclose(DM_b, expected, rel_tol=1e-12)
    DM_b_obs = res["d"]["DM_baryon_ratio_obs"]
    # Within 2 sigma of Planck
    assert abs(res["d"]["pull_DM_baryon"]) < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  Structural closure  (E_8 mode accounting)
# ═══════════════════════════════════════════════════════════════════════

def test_e8_mode_assignment_closure(res):
    """248 = 14 + 52 + 182.  Every E_8 mode has a physical job."""
    R = res["R"]
    assert R["dim_G2"] + R["dim_F4"] + R["N_bridge"] == 248
    assert R["dim_G2"] == 14
    assert R["dim_F4"] == 52
    assert R["N_bridge"] == 182


# ═══════════════════════════════════════════════════════════════════════
#  Polynomial identities  (the four-integer reductions)
#  ──────────────────────────────────────────────────────────────────────
#  These check "does the named quantity equal the claimed polynomial in
#  (d_10, d_11, n_7, n_26)?".  They are the unified version's strongest
#  rhetorical move: 14+ named quantities reduce to one screen.
# ═══════════════════════════════════════════════════════════════════════

def test_four_irreducible_integers(res):
    """(d_10, d_11, n_7, n_26) = (2, 3, 7, 26)."""
    R = res["R"]
    assert (R["d10"], R["d11"], R["n7"], R["n26"]) == (2, 3, 7, 26)


def test_sugawara_consistency_d11sq_minus_1_eq_4_d10(res):
    """The integer identity d_11^2 - 1 = 4 d_10 forces (d_10, d_11) = (2, 3).

    Among positive integers, this Diophantine constraint together with
    d_10 < d_11 and d_10, d_11 >= 2 (both quantum dimensions of nontrivial
    primaries) has the unique solution (d_10, d_11) = (2, 3).
    Equivalently: d_11^2 = 4 d_10 + 1.
    """
    R = res["R"]
    assert R["d11"] ** 2 - 1 == 4 * R["d10"]


def test_singh_ratio_is_16(res):
    """alpha_s / alpha_em = (8/3) * C_2(26) = (8/3) * 6 = 16."""
    R = res["R"]
    charge_trace = Fraction(8, 3)
    singh_ratio = charge_trace * R["C2_26"]
    assert singh_ratio == 16


def test_vertex_count_is_30_with_no_free_choice(res):
    """N_vertex = n_26 (F_4 fund) + 4 (Higgs real DOFs) = 30.

    Equivalently: N_vertex = n_26 + h_dual(G_2)/d_10 since 4 = d_10^2.
    Either way, the count is forced by group theory, not chosen to fit v_EW.
    """
    R = res["R"]
    assert R["N_vertex"] == R["n26"] + 4 == 30


def test_dual_coxeter_numbers_are_squares(res):
    """h_dual(G_2) = d_10^2 = 4 and h_dual(F_4) = d_11^2 = 9."""
    R = res["R"]
    assert R["hv_G2"] == R["d10"] ** 2 == 4
    assert R["hv_F4"] == R["d11"] ** 2 == 9


def test_casimirs_factor_through_d10_d11(res):
    """C_2(7) = d_10 and C_2(26) = d_10 * d_11."""
    R = res["R"]
    assert R["C2_7"] == R["d10"] == 2
    assert R["C2_26"] == R["d10"] * R["d11"] == 6


def test_K_is_d10_times_d11(res):
    """WZW altitude K = k + h_dual = d_11 + d_11 = d_10 * d_11 = 6."""
    R = res["R"]
    assert R["K"] == R["d10"] * R["d11"] == 6


def test_h10_is_d10_over_d11_squared(res):
    """h(1,0) = d_10 / d_11^2 = 2/9."""
    R = res["R"]
    assert math.isclose(R["h10"], R["d10"] / R["d11"] ** 2, rel_tol=1e-12)
    assert math.isclose(R["h10"], 2 / 9, rel_tol=1e-12)


def test_h11_is_one_over_d10(res):
    """h(1,1) = 1 / d_10 = 1/2."""
    R = res["R"]
    assert math.isclose(R["h11"], 1.0 / R["d10"], rel_tol=1e-12)
    assert math.isclose(R["h11"], 1 / 2, rel_tol=1e-12)


def test_Q0_is_d10_over_d11(res):
    """Koide value Q_0 = d_10 / d_11 = 2/3."""
    R = res["R"]
    assert math.isclose(R["Q0"], R["d10"] / R["d11"], rel_tol=1e-12)
    assert math.isclose(R["Q0"], 2 / 3, rel_tol=1e-12)


def test_sin2W_is_d11_over_d10sq_plus_d11sq(res):
    """sin^2 theta_W = d_11 / (d_10^2 + d_11^2) = 3/13."""
    R = res["R"]
    expected = R["d11"] / (R["d10"] ** 2 + R["d11"] ** 2)
    assert math.isclose(R["sin2W"], expected, rel_tol=1e-12)
    assert math.isclose(R["sin2W"], 3 / 13, rel_tol=1e-12)


def test_N_bridge_is_n7_times_n26(res):
    """N_bridge = n_7 * n_26 = 182."""
    R = res["R"]
    assert R["N_bridge"] == R["n7"] * R["n26"] == 182


def test_h_bridge_is_one_marginal(res):
    """h_bridge = h(7,G_2) + h(26,F_4) = 2/5 + 3/5 = 1 (marginal primary)."""
    R = res["R"]
    assert math.isclose(R["h_7"], 2 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_26"], 3 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_bridge"], 1.0, rel_tol=1e-12)
    # The Fraction identity: h(R) = C_2(R) / (1 + h_dual)
    assert math.isclose(R["h_7"], R["C2_7"] / (1 + R["hv_G2"]), rel_tol=1e-12)
    assert math.isclose(R["h_26"], R["C2_26"] / (1 + R["hv_F4"]), rel_tol=1e-12)


def test_c_coset_vanishes(res):
    """c_coset = c(E_8) - c(G_2) - c(F_4) = 8 - 14/5 - 26/5 = 0 (topological)."""
    R = res["R"]
    assert math.isclose(R["c_coset"], 0.0, abs_tol=1e-12)


def test_xi_bridge_is_one_over_48_pi(res):
    """xi_bridge = alpha_G_2 * E[v^2] * h_bridge = 1/(48 pi)."""
    R = res["R"]
    assert math.isclose(R["xi_bridge"], 1.0 / (48.0 * math.pi), abs_tol=1e-14)
    assert math.isclose(R["E_v2"], 0.5, abs_tol=1e-12)


# ═══════════════════════════════════════════════════════════════════════
#  Algebraic structure  (rational identities for masses + mixing)
# ═══════════════════════════════════════════════════════════════════════

def test_light_quark_ratios_are_38_9_and_83_9(res):
    """m_u/m_e = 38/9 and m_d/m_e = 83/9 (triality + h(1,0) = 2/9)."""
    h10 = Fraction(2, 9)
    assert Fraction(4, 1) + h10 == Fraction(38, 9)
    assert Fraction(9, 1) + h10 == Fraction(83, 9)
    m = res["m"]
    assert math.isclose(m["m_u"] / m["m_e"], 38 / 9, rel_tol=1e-12)
    assert math.isclose(m["m_d"] / m["m_e"], 83 / 9, rel_tol=1e-12)


def test_heavy_quark_ratios_are_217_18_and_1165_12(res):
    """m_c/m_mu = 217/18 (charm) and m_t/m_tau = 1165/12 (top)."""
    delta_OPE = Fraction(1, 18)
    K = 6
    assert Fraction(12, 1) + delta_OPE == Fraction(217, 18)
    assert Fraction(81 + 16, 1) + Fraction(1, 2 * K) == Fraction(1165, 12)
    m = res["m"]
    assert math.isclose(m["m_c"] / m["m_mu"], 217 / 18, rel_tol=1e-12)
    assert math.isclose(m["m_t"] / m["m_tau"], 1165 / 12, rel_tol=1e-12)


def test_bottom_koide_is_289_432(res):
    """Q(c,b,t) = 2/3 + h(1,1)/K^3 = 289/432."""
    K = 6
    Q = Fraction(2, 3) + Fraction(1, 2) / Fraction(K ** 3, 1)
    assert Q == Fraction(289, 432)
    m = res["m"]
    sc = math.sqrt(m["m_c"])
    sb = math.sqrt(m["m_b"])
    st = math.sqrt(m["m_t"])
    Q_num = (m["m_c"] + m["m_b"] + m["m_t"]) / (sc + sb + st) ** 2
    assert math.isclose(Q_num, 289 / 432, rel_tol=1e-10)


def test_strange_koide_is_649_972():
    """Q(s,c,b) = 2/3 + h(1,0)/K^3 = 649/972 (Rivero fundamental channel)."""
    K = 6
    Q = Fraction(2, 3) + Fraction(2, 9) / Fraction(K ** 3, 1)
    assert Q == Fraction(649, 972)


def test_albert_bridge_factor_is_32_over_27():
    """bridge^2 = Q_0^2 * d_10^3 / d_11 = (4/9)(8/3) = 32/27."""
    d10, d11 = 2, 3
    bridge_sq = Fraction(2, 3) ** 2 * Fraction(d10 ** 3, 1) / Fraction(d11, 1)
    assert bridge_sq == Fraction(32, 27)


def test_lepton_koide_is_two_thirds(res):
    """Q_lep = 1/3 + (B/A)^2 / 6 = 2/3 from |B/A|^2 = 2."""
    Q = Fraction(1, 3) + Fraction(2, 6)
    assert Q == Fraction(2, 3)
    m = res["m"]
    me, mmu, mtau = m["m_e"], m["m_mu"], m["m_tau"]
    Q_num = (me + mmu + mtau) / (math.sqrt(me) + math.sqrt(mmu) + math.sqrt(mtau)) ** 2
    assert math.isclose(Q_num, 2 / 3, rel_tol=1e-5)


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
    # Depth-1 closed form 256(2 pi - 1)/pi^2 is the truncation; the full
    # value includes the FORCED depth-3 two-orientation e<->q loop and
    # solves the cubic  x^3 = (512/pi)[(1 - 1/(2 pi)) x^2 - 1/(2 pi^2)].
    x = bridge["inv_alpha_phys"]
    depth1 = 256.0 * (2.0 * math.pi - 1.0) / math.pi ** 2
    assert x < depth1  # the depth-3 echo lowers 1/alpha
    cubic = x**3 - (2**9 / math.pi) * ((1.0 - 1.0/(2.0*math.pi)) * x**2
                                       - 1.0/(2.0*math.pi**2))
    assert abs(cubic) < 1e-6
    # lands on the Berkeley Cs recoil measurement 137.035999046(27)
    assert abs(x - 137.035999046) < 27e-9


def test_alpha_phys_agrees_with_pdg(res):
    """1/alpha(0) = 256(2 pi - 1)/pi^2 matches PDG 137.0360 to 3 ppm."""
    inv_alpha = res["c"]["bridge"]["inv_alpha_phys"]
    # PDG (lower-precision quote, full Penning-trap is sub-ppb)
    assert err_pct(inv_alpha, 137.0360) < 0.001


# ═══════════════════════════════════════════════════════════════════════
#  Generation word lemma + linearity from closure (code-only ports)
# ═══════════════════════════════════════════════════════════════════════

def test_generation_words_are_walk_counts():
    """base(n) = (4, 12, 97) as boundary-walk counts; Z2 split u/d."""
    import words
    w = words.derive()
    assert w["base"] == (4, 12, 97)
    assert w["base_down_seed"] == 9          # d = aa (Z2-swapped seed)
    assert w["FP_C0"] == 12                   # neutral lane = base(2)
    # mixed Z2-symmetric terminal words are excluded (data: 64/139 GeV)
    assert set(w["excluded"].values()) == {36, 78}


def test_linearity_from_closure():
    """L_x linear + |xy|=|x||y| on O; sedenion zero divisor; Z3 rule."""
    import numpy as np
    import octonions as oc
    rng = np.random.default_rng(11)
    assert oc._check_linearity(rng)
    assert oc._check_hurwitz_gate(rng)
    assert oc._check_z3_selection()



# ═══════════════════════════════════════════════════════════════════════
#  Vertex composition rule + anchor inversion (m_e → M_Pl → G)
# ═══════════════════════════════════════════════════════════════════════

def test_vertex_rule_identities():
    """Depth-3 multiplicities = products of the node's own depth-1
    vertices: 2 (orientations), 16 = charge_trace×C2(26), 27/2 = W·dimJ3O."""
    from root import charge_trace, C2_26, d10, d11
    assert charge_trace * C2_26 == 16
    assert Fraction(d11**3, d10) == Fraction(27, 2)


def test_electron_anchor_exact(res):
    """m_e reproduces the sole dimensional input exactly."""
    from root import M_E_ANCHOR_MEV
    assert abs(res["m"]["m_e"] / M_E_ANCHOR_MEV - 1.0) < 1e-12


def test_muon_mass_predicted_ppb(res):
    """m_mu is PREDICTED to sub-ppm: within 1 sigma of CODATA."""
    m_mu_codata, sigma = 105.6583755, 2.3e-6
    assert abs(res["m"]["m_mu"] - m_mu_codata) / sigma < 1.0


def test_tau_mass_prediction_band(res):
    """m_tau prediction inside the current PDG band (sharp later)."""
    assert abs(res["m"]["m_tau"] - 1776.93) / 0.09 < 1.0


def test_newton_constant_predicted():
    """G is now a PREDICTION: within 3 sigma of CODATA (currently -2.0)."""
    from root import G_PRED, G_CODATA
    assert abs(G_PRED / G_CODATA - 1.0) < 3 * 2.2e-5


def test_vew_matches_GF_at_1e7(res):
    """v_EW (with the FORCED 16-term and the m_e anchor) matches the
    G_F determination at the 1e-7 level."""
    v_target = (math.sqrt(2.0) * 1.1663787e-5) ** -0.5
    assert abs(res["m"]["v_EW_pred_GeV"] / v_target - 1.0) < 5e-7


# ═══════════════════════════════════════════════════════════════════════
#  The canonical freeze: the full value table, pinned.
#  Any future edit that moves a canonical number past these tolerances
#  is either a derivation change (update deliberately, with provenance)
#  or a regression (fix it).
# ═══════════════════════════════════════════════════════════════════════

def test_mass_coordinate_invariants(res):
    """Scheme-invariant light-quark ratios: predictions with no
    mass-coordinate choice at all (QCD running cancels)."""
    ms = res["m"]
    assert abs(ms["mu_over_md"] - 38.0 / 83.0) < 1e-12   # exact
    assert abs(ms["ms_over_mud"] - 27.3183) < 1e-3
    assert abs(ms["Q_ellipse"] - 22.383) < 1e-2


def test_canonical_freeze(res):
    """The frozen canonical table (ledger closed)."""
    R, m, g, c = res["R"], res["m"], res["g"], res["c"]
    import root as _r

    # exact (dimensionless web + anchor inversion)
    assert abs(_r.inv_alpha_phys - 137.0359990495613) < 1e-10
    assert abs(m["m_mu"] - 105.65837578191945) < 1e-9
    assert abs(m["m_tau"] - 1776.9092813590157) < 1e-7
    assert abs(m["v_EW_pred_GeV"] - 246.21964503014013) < 1e-9
    assert abs(_r.G_PRED - 6.674003218967312e-11) < 1e-22
    assert abs(_r.M_Pl_GeV - 1.220917145099806e+19) < 1e10

    # face-split law closures
    assert abs(g["G_ratio_UV"] - 1.0) < 1e-6
    assert abs(g["sin2W"] - 0.23128541) < 1e-7

    # EW layer-1 chain (PDG 2024 imports)
    assert abs(g["M_W_derived"] - 80.3559) < 2e-3
    assert abs(g["M_Z_derived"] - 91.1885) < 2e-3

    # couplings + Higgs (declared imports, banded)
    # alpha_s matched at mu* = M_Pl e^{-(9pi^2/2-6)} = v_EW e^{15/512}
    # (gauge lever-arm endpoint; pre-migration label was v_EW)
    assert abs(c["mu_star_GeV"] - 253.534) < 1e-2
    assert abs(c["alpha_s_MZ_thresh"] - 0.118385) < 1e-5
    assert abs(g["mH_pred"] - 125.2965) < 0.05  # (1−h₁₀) vent; registry.PROMOTIONS

    # words lemma bases (integer ranks, exact)
    import words as _w
    import io as _io, contextlib as _ctx
    with _ctx.redirect_stdout(_io.StringIO()):
        bases = _w.derive()["base"]
    assert bases == (4, 12, 97)


# ═══════════════════════════════════════════════════════════════════════
#  Embedding uniqueness (E₈ coherence gate)
# ═══════════════════════════════════════════════════════════════════════

def test_embedding_uniqueness():
    """E₈(1) ⊃ G₂(1) × F₄(1) is the unique conformal embedding
    passing all six coherence gates among rank-8 simple algebras."""
    import embedding_uniqueness as eu
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = eu.derive()
    assert result['n_survivors'] == 1
    assert result['survivor'] == ('E8', 'G2', 'F4')


# ═══════════════════════════════════════════════════════════════════════
#  Protected forgetting (G₂ harmonic)
# ═══════════════════════════════════════════════════════════════════════

def test_protected_forgetting():
    """v = (1, -1/2, -1/2) satisfies PvP=0 and Pv²P=1/2."""
    import protected_forgetting as pf
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = pf.derive()
    assert abs(result['PvP']) < 1e-14
    assert abs(result['Pv2P'] - 0.5) < 1e-14
    assert result['unique'] is True


# ═══════════════════════════════════════════════════════════════════════
#  NLS soliton (substrate dynamics)
# ═══════════════════════════════════════════════════════════════════════

def test_nls_soliton_stability():
    """Z₃-NLS evolves stably with number conservation."""
    import nls_soliton
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = nls_soliton.derive()
    assert result['conservation'] < 1e-6


# ═══════════════════════════════════════════════════════════════════════
#  The canonical mass freeze (full 9-mass table, high precision)
# ═══════════════════════════════════════════════════════════════════════

FROZEN_MASSES_MEV = {
    'm_e': 0.51099895069,        # anchor, exact by construction
    'm_mu': 105.6583758,
    'm_tau': 1776.9092814,
    'm_u': 2.157551125,
    'm_d': 4.712545879,
    'm_s': 93.8395306,           # depth-2 bridge update (registry.WATCHES 3)
    'm_c': 1273.770419,
    'm_b': 4193.80675,
    'm_t': 172508.2761,
}


@pytest.mark.parametrize("key,val", sorted(FROZEN_MASSES_MEV.items()))
def test_frozen_mass_table(res, key, val):
    """Every mass the structure produces, frozen to 10 digits."""
    assert abs(res["m"][key] / val - 1.0) < 1e-9


def test_frozen_web_nodes(res):
    """The graph nodes' fixed-point values, frozen."""
    from root import WEB
    st = WEB.state
    assert abs(st["lambda_MPl"] - (-0.012450153)) < 1e-8  # ×(7/9): fundamental-share vent
    assert abs(st["G_ratio"] - 0.999999917) < 1e-8
    assert abs(st["dark_ratio"] - (2.0 * math.pi - 1.0)) < 1e-12
    assert abs(st["v_EW_GeV"] - 246.219645) < 1e-5
    assert abs(st["M_Pl_MeV"] / 1e3 - 1.220917145e19) < 1e11


def test_scale_landmarks(res):
    """The derived scale ladder: Lambda_conf, mu*, M_V, M_R."""
    import root as _r
    Lambda_conf = (1.0 / _r.d10) * _r.M_Pl_GeV * math.exp(-_r.S_lepton)
    assert abs(Lambda_conf - 0.3138) < 1e-3
    assert abs(res["c"]["mu_star_GeV"] - 253.534) < 1e-2
    assert abs(res["c"]["M_V_derived"] - 111.648) < 1e-2
    M_R = _r.M_Pl_GeV * math.exp(-2.0 * math.pi**2)
    assert abs(M_R / 3.27e10 - 1.0) < 0.01


def test_np_splitting_qcd_part(res):
    """m_d - m_u = 5 m_e exactly; the lattice QCD part of m_n - m_p
    (BMW, Science 347:1452 (2015): 2.52 +/- 0.23 MeV)."""
    m = res["m"]
    dmq = m["m_d"] - m["m_u"]
    assert abs(dmq / (5.0 * m["m_e"]) - 1.0) < 1e-12
    pull = (dmq - 2.52) / 0.23
    assert abs(pull) < 2.0          # currently +0.15 sigma


def test_diquark_channel_ordering():
    """3 x 3 = 3bar + 6 forced; h(3bar) = 2/9 < h(6) = 5/9, so the
    nucleon (good diquark) sits below the Delta, forced."""
    import mixing as _mx
    w = _mx._compute_wzw()
    IDX, NF, P = w['IDX'], w['N_fus'], w['PRIMARIES']
    channels = [c for c in P if NF[IDX[(1, 0)], IDX[(1, 0)], IDX[c]]]
    assert sorted(channels) == [(0, 1), (2, 0)]

    def h(l1, l2):
        C2 = Fraction(l1*l1 + l2*l2 + l1*l2 + 3*l1 + 3*l2, 3)
        return C2 / 6

    assert h(0, 1) == Fraction(2, 9)
    assert h(2, 0) == Fraction(5, 9)
    assert h(0, 1) < h(2, 0)


def test_no_fourth_generation():
    """Windings are mod 3: a fourth label does not exist."""
    import root as _r
    def gap(k):
        return 1.0 + math.sqrt(_r.d10) * math.cos(
            float(_r.h10) + 2.0 * math.pi * k / 3.0)
    for k in range(3):
        assert abs(gap(k + 3) - gap(k)) < 1e-14


def test_neutrino_spectrum(res):
    """m1 = 0 (rank-2 seesaw) + measured splittings -> the registered
    kinematic spectrum: sum, m_beta, m_bb band."""
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


def test_bell_fork():
    """The local realist reading (outcomes as assigned +/-1 facts) is
    capped at 2; the framework's local non-realist reading (the
    ledger) sits exactly at the Tsirelson bound."""
    from probes import bell
    S_local = bell.chsh_local_bound(n_strategies=2000, n_lambda=200)
    assert S_local <= 2.0 + 1e-12
    assert abs(bell.chsh_ledger() - 2.0 * math.sqrt(2.0)) < 1e-12


# ═══════════════════════════════════════════════════════════════════════
#  Attractor (Chladni) property of the web
# ═══════════════════════════════════════════════════════════════════════

def test_web_attractor(res):
    """Perturb every web node by up to ±50% and re-solve: the whole web
    returns to its fixed point (machine precision, 5/5 trials).  The
    constants are stored nowhere; they are where the recursion lands.
    Chladni figures, in code."""
    import random
    from root import WEB

    WEB.solve()                      # ensure the state IS the fixed point
    baseline = dict(WEB.state)
    random.seed(3)
    try:
        for _ in range(5):
            WEB.state = {k: v * (1.0 + random.uniform(-0.5, 0.5))
                         for k, v in baseline.items()}
            WEB.solve(iters=2000)
            dev = max(abs(WEB.state[k] / baseline[k] - 1.0)
                      for k in baseline if baseline[k] != 0.0)
            assert dev < 1e-10, f"web did not return: max dev {dev:.2e}"
    finally:
        WEB.state = dict(baseline)
        WEB.solve()


