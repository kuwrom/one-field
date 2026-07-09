"""Scorecard tests for the interference formulation.

Every numerical prediction and every polynomial identity is
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
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Mass sector  (1 anchor + 8 predictions)
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("name,tol_pct", [
    ("e",   0.01),
    ("mu",  0.01),
    ("tau", 0.01),
])
def test_lepton_mass(res, name, tol_pct):
    """Completed charged-lepton Web readings from the Z_3 circulant."""
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


def test_full_record_strong_cp_closure(res):
    """The production chain exposes the complete-vacuum character,
    not separately basis-dependent bare-theta and mass phases."""
    sc = res["c"]["strong_cp"]
    assert res["c"]["bar_theta"] == 0.0
    assert sc["physical_angle"] == "bar_theta"
    assert sc["phase_rule_status"] == "CONSTITUTIVE ONTOLOGY"
    assert sc["vacuum_phase_automorphisms"] == (1.0 + 0j,)
    assert all(value == 1.0 + 0j
               for value in sc["winding_character_samples"].values())
    assert not sc["theta_QCD_separately_defined"]
    assert not sc["mass_determinant_phase_separately_defined"]


# ═══════════════════════════════════════════════════════════════════════
#  CKM sector  (Wolfenstein + UT angles, summarised)
# ═══════════════════════════════════════════════════════════════════════

def test_ckm_diagonal_pull_score(res):
    """CKM correlated coordinates have a small descriptive pull score."""
    n = len(res["mx"]["ckm_obs"])
    score_per_coordinate = res["mx"]["diagonal_pull_score_ckm"] / n
    assert score_per_coordinate < 1.0


def test_ckm_max_pull(res):
    """No individual CKM observable deviates by more than 2 sigma."""
    assert res["mx"]["max_pull_ckm"] < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  PMNS sector (3 mixing angles + convention-dependent U_e3 phase)
# ═══════════════════════════════════════════════════════════════════════

def test_pmns_diagonal_pull_score(res):
    """PMNS marginal pulls have a small descriptive diagonal score."""
    score_per_angle = res["mx"]["diagonal_pull_score_pmns"] / 3
    assert score_per_angle < 1.0


def test_pmns_generation_block_is_structural(res):
    """The projected bulk-fusion blocks have the stated TBM pattern."""
    assert res["mx"]["generation_orbits"] == (
        ((0, 1), (1, 0)), ((0, 2), (2, 0)), ((1, 2), (2, 1)))
    assert res["mx"]["democratic_eigenvalues"] == pytest.approx(
        (3.0, 0.0, 0.0), abs=1e-12)
    assert res["mx"]["fundamental_eigenvalues"] == pytest.approx(
        (-0.5, 1.0, 0.5), abs=1e-12)
    tbm_abs = np.array([
        [2/math.sqrt(6), 1/math.sqrt(3), 0],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
    ])
    assert np.allclose(np.abs(res["mx"]["U_TBM"]), tbm_abs)


def test_pmns_fold_is_not_mislabeled_as_the_zc_nimrep(res):
    """Z_C has two boundaries; the three-orbit projection is not closed."""
    assert res["mx"]["zc_boundary_states"] == 2
    assert np.array_equal(res["mx"]["zc_n10"], np.ones((2, 2), dtype=int))
    assert np.array_equal(
        res["mx"]["zc_n11"], np.array([[1, 2], [2, 1]], dtype=int))
    assert res["mx"]["zc_fusion_residual"] == 0.0
    assert res["mx"]["fundamental_generation_leakage"] > 1e-12
    assert res["mx"]["fold_fusion_residual"] > 1e-12


def test_ckm_d6_boundary_data_are_computed(res):
    """The pair fraction and braid magnitude come from the graph/cells."""
    assert res["mx"]["d6_boundary_pairs"] == ((0, 5), (1, 2), (3, 4))
    expected = np.array([[1.0, 2.0], [2.0, 1.0]])
    assert all(np.array_equal(block, expected)
               for block in res["mx"]["d6_adjoint_blocks"])
    assert res["mx"]["A"]**2 == pytest.approx(2/3, abs=1e-15)
    cell_hi, cell_lo = res["mx"]["d6_cell_magnitudes"]
    assert cell_hi**2 == pytest.approx((math.sqrt(3) + 1)/2, abs=1e-15)
    assert cell_lo**2 == pytest.approx((math.sqrt(3) - 1)/2, abs=1e-15)
    assert res["mx"]["d6_connection_magnitude"] == pytest.approx(
        1/math.sqrt(2), abs=1e-15)
    assert res["mx"]["d6_braid_residual"] < 1e-12
    assert sorted(res["mx"]["d6_depth4_carrier_dimensions"]) == [2, 3, 3]


def test_pmns_ue3_phase_output(res):
    """The identified PMNS map returns a finite phase-convention output."""
    ue3_phase = res["mx"]["ue3_phase_deg"]
    assert 0 < ue3_phase < 180


# ═══════════════════════════════════════════════════════════════════════
#  Higgs running coordinate
# ═══════════════════════════════════════════════════════════════════════

def test_higgs_transport_coordinate(res):
    """The candidate boundary is transported in one MS-bar scheme."""
    g = res["g"]
    assert g["lambda_mt_transport"] == pytest.approx(0.1265078, abs=2e-7)
    assert g["higgs_transport_coordinate_GeV"] == pytest.approx(
        123.8501, abs=2e-3)


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
#  Gravity sector  (Newton constant from the (7,26) bridge heat kernel)
# ═══════════════════════════════════════════════════════════════════════

def test_newton_internal_closure_UV(res):
    """Freeze the gravity normalisation's internal UV closure; not an independent G test."""
    assert abs(res["g"]["G_ratio_UV"] - 1.0) < 1e-6


def test_newton_internal_closure_broken(res):
    """Freeze the broken-phase diagnostic alongside the gravity normalisation's UV value."""
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
#  Electroweak (derived)  (sin^2 theta_W, M_Z, M_W)
# ═══════════════════════════════════════════════════════════════════════

def test_sin2_theta_W(res):
    """sin^2 theta_W = 3/13 + Q₀²·α/(2π), depth-1 h₇ + depth-2 h₇/d₁₁²."""
    sin2W = res["g"]["sin2W"]
    # canonical: base 3/13 + depth-1 h₇·α/(2π) + depth-2 (h₇/d₁₁²)·α/(2π)
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
    from interference import wiring

    DM_b = res["d"]["DM_baryon_ratio"]
    expected = 2.0 * math.pi - 1.0
    assert math.isclose(DM_b, expected, rel_tol=1e-12)
    assert res["d"]["f_baryonic"] == wiring.normalized_current_response(1.0)
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


def test_vertex_current_trace_is_30(res):
    """The E8 positive-root current index supplies the vertex factor."""
    R = res["R"]
    assert R["vertex_charge_trace"] == R["N_vertex"] == R["hv_E8"] == 30


#  Constants root.py DEFINES by exactly the expression a per-constant
#  test would restate (hv_G2 = d10**2, K = d10*d11, h10 =
#  Fraction(d10, d11**2), ...).  Re-deriving the expression tests
#  nothing; what carries regression value is the VALUE, so an edit to
#  any of root.py's definitions fails here.  The identities that cross
#  two INDEPENDENT routes are not in this table and keep their own
#  tests below: Q0 (defined as 1/3 + d10/6, asserted equal to
#  d10/d11), h_bridge (a sum asserted exactly one), c_coset (a
#  difference asserted exactly zero), xi_bridge (a product of three
#  separately derived factors landing on 1/(48 pi)).
_ROOT_CONSTANTS = [
    ("hv_G2", 4), ("hv_F4", 9),         # dual Coxeter numbers
    ("C2_7", 2), ("C2_26", 6),          # Casimirs
    ("K", 6),                           # WZW altitude
    ("h10", 2 / 9), ("h11", 1 / 2),     # conformal weights
    ("sin2W", 3 / 13),                  # zero-edge Weinberg base
    ("N_bridge", 182),                  # bridge mode count
]


def test_root_constants_frozen(res):
    """The defined constants of the algebraic root, frozen by value."""
    R = res["R"]
    for key, expected in _ROOT_CONSTANTS:
        assert math.isclose(float(R[key]), expected, rel_tol=1e-12), key


def test_Q0_is_d10_over_d11(res):
    """Koide value Q_0 = d_10 / d_11 = 2/3."""
    R = res["R"]
    assert math.isclose(R["Q0"], R["d10"] / R["d11"], rel_tol=1e-12)
    assert math.isclose(R["Q0"], 2 / 3, rel_tol=1e-12)


def test_h_bridge_is_one(res):
    """h_bridge = h(7,G_2) + h(26,F_4) = 2/5 + 3/5 = 1 (dimension one)."""
    R = res["R"]
    assert math.isclose(R["h_7"], 2 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_26"], 3 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_bridge"], 1.0, rel_tol=1e-12)
    # The Fraction identity: h(R) = C_2(R) / (1 + h_dual)
    assert math.isclose(R["h_7"], R["C2_7"] / (1 + R["hv_G2"]), rel_tol=1e-12)
    assert math.isclose(R["h_26"], R["C2_26"] / (1 + R["hv_F4"]), rel_tol=1e-12)


def test_c_coset_vanishes(res):
    """c_coset = c(E_8) - c(G_2) - c(F_4) = 8 - 14/5 - 26/5 = 0."""
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
    # Depth-1 closed form 256(2 pi - 1)/pi^2 is the truncation. The full
    # value includes the FORCED depth-4 reading of the e<->q carrier and
    # solves the cubic  x^3 = (512/pi)[(1 - 1/(2 pi)) x^2 - 1/(2 pi^2)].
    x = bridge["inv_alpha_phys"]
    depth1 = 256.0 * (2.0 * math.pi - 1.0) / math.pi ** 2
    assert x < depth1  # the depth-4 response lowers 1/alpha
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
    """base(n) = (4, 12, 97) as boundary-walk counts. Z2 split u/d."""
    import words
    w = words.derive()
    assert w["base"] == (4, 12, 97)
    assert w["base_down_seed"] == 9          # d = aa (Z2-swapped seed)
    assert w["FP_C0"] == 12                   # neutral lane = base(2)
    # mixed Z2-symmetric terminal words are excluded (data: 64/139 GeV)
    assert set(w["excluded"].values()) == {36, 78}


def test_linearity_from_closure():
    """L_x linear + |xy|=|x||y| on O. Sedenion zero divisor. Z3 rule."""
    import numpy as np
    import octonions as oc
    rng = np.random.default_rng(11)
    assert oc._check_linearity(rng)
    assert oc._check_hurwitz_gate()
    assert oc._check_z3_selection()


def test_sedenion_zero_divisor_census():
    """The '168 pairs' remark is computed, not quoted: full enumeration
    of (e_i ± e_j)(e_k ± e_l) = 0 on the unit-coefficient sedenion
    lattice gives 336 ordered = 168 unordered zero-divisor pairs."""
    import octonions as oc
    assert oc._check_hurwitz_gate(census=True)



# ═══════════════════════════════════════════════════════════════════════
#  Vertex composition rule + anchor inversion (m_e → M_Pl → G)
# ═══════════════════════════════════════════════════════════════════════

def test_vertex_rule_identities():
    """Depth-4 responses = products of the node's own depth-1
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
    """G is a PREDICTION: within 3 sigma of CODATA (-2.0)."""
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
    m, g, c = res["m"], res["g"], res["c"]
    import root as _r

    # exact (dimensionless web + anchor inversion)
    assert abs(_r.inv_alpha_phys - 137.0359990495613) < 1e-10
    assert abs(m["m_mu"] - 105.65837578191945) < 1e-9
    assert abs(m["m_tau"] - 1776.9092813590157) < 1e-7
    assert abs(m["v_EW_pred_GeV"] - 246.21964503014013) < 1e-9
    assert abs(_r.G_PRED - 6.674003218967312e-11) < 1e-22
    assert abs(_r.M_Pl_GeV - 1.220917145099806e+19) < 1e10

    # apportionment law closures
    assert abs(g["G_ratio_UV"] - 1.0) < 1e-6
    assert abs(g["sin2W"] - 0.23128541) < 1e-7

    # EW transport chain (PDG 2024 imports)
    assert abs(g["M_W_derived"] - 80.3559) < 2e-3
    assert abs(g["M_Z_derived"] - 91.1885) < 2e-3

    # couplings + Higgs (declared imports, banded)
    # alpha_s matched at mu* = M_Pl e^{-(9pi^2/2-6)} = v_EW e^{15/512}
    # (gauge lever-arm endpoint. Pre-migration label was v_EW)
    assert abs(c["mu_star_GeV"] - 253.534) < 1e-2
    assert abs(c["alpha_s_MZ_thresh"] - 0.118385) < 1e-5
    assert abs(g["lambda_mt_transport"] - 0.1265078) < 2e-7

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
    """E₈(1) ⊃ G₂(1) × F₄(1) is the unique pair passing all six
    coherence gates among the rank-8 candidates, computed from the
    scan, not asserted as literals, with the E₇ × A₁ contender
    (passes (a)-(e), fails only trivial centres) frozen alongside:
    gates (e) and (f) are framework closure requirements, and the scan
    keeps their decisive role visible."""
    import embedding_uniqueness as eu
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = eu.derive()
    assert result['n_survivors'] == 1
    assert result['survivor'] == ('E8', 'G2', 'F4')
    nm = {(p, frozenset((g, f))): failed
          for p, g, f, failed in result['near_misses']}
    assert nm == {('E8', frozenset(('E7', 'A1'))): '(f) trivial centres'}


# ═══════════════════════════════════════════════════════════════════════
#  Protected forgetting (G₂ harmonic)
# ═══════════════════════════════════════════════════════════════════════

def test_protected_forgetting():
    """The canonical Z₃ harmonic (B, θ) = (√2, 2/9): PvP = 0, total
    second moment Pv²P = 1, per-mode power E[v²] = 1/2, Koide Q = 2/3.

    The retired harmonic (1, -1/2, -1/2) is (B, θ) = (1, 0), which gives
    Q = 1/2 with two masses degenerate; θ = 0 is degenerate for every
    amplitude.  Both facts are computed, not asserted.
    """
    import protected_forgetting as pf
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = pf.derive()
    assert abs(result['PvP']) < 1e-14
    # total vs per-mode: Pv²P = 2·E[v²], and the web consumes E[v²]
    assert abs(result['Pv2P_total'] - 1.0) < 1e-12
    assert abs(result['E_v2_mode'] - 0.5) < 1e-12
    assert abs(result['Pv2P_total'] - 2 * result['E_v2_mode']) < 1e-12
    # Koide value, and the amplitude it pins
    assert abs(result['Q'] - 2.0 / 3.0) < 1e-12
    assert result['amplitude_pinned_by_koide'] is True
    # the retired harmonic really is degenerate at Q = 1/2
    assert abs(result['retired_Q'] - 0.5) < 1e-12
    assert result['theta0_degenerate'] is True
    # the short-root direction gives three distinct
    assert result['short_root_three_distinct'] is True
    # canonical harmonic stays on the positive-gap branch
    assert result['gap_min'] > 0


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
    'm_s': 93.8395306,           # depth-2 bridge update
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
    """The scale ladder, with M_V the connector's per-mode reading."""
    import root as _r
    Lambda_conf = (1.0 / _r.d10) * _r.M_Pl_GeV * math.exp(-_r.S_lepton)
    assert abs(Lambda_conf - 0.3138) < 1e-3
    assert abs(res["c"]["mu_star_GeV"] - 253.534) < 1e-2
    M_V = res["c"]["M_V_GeV"]
    assert abs(M_V - 111.648) < 1e-2
    g_G2 = math.sqrt(4.0 * math.pi * res["c"]["alpha_s_mu_star"])
    n_broken = 2 * _r.n7 - (_r.d11**2 - 1)
    v_break = res["m"]["v_EW_pred_GeV"] / math.sqrt(_r.d10)
    assert M_V == pytest.approx(
        g_G2 * v_break * math.sqrt(_r.C2_7 / n_broken), abs=1e-12)
    M_R = _r.M_Pl_GeV * math.exp(-2.0 * math.pi**2)
    assert abs(M_R / 3.27e10 - 1.0) < 0.01


def test_np_splitting_qcd_part(res):
    """m_d - m_u = 5 m_e exactly. The lattice QCD part of m_n - m_p
    (BMW, Science 347:1452 (2015): 2.52 +/- 0.23 MeV)."""
    m = res["m"]
    dmq = m["m_d"] - m["m_u"]
    assert abs(dmq / (5.0 * m["m_e"]) - 1.0) < 1e-12
    pull = (dmq - 2.52) / 0.23
    assert abs(pull) < 2.0          # currently +0.15 sigma


def test_diquark_channel_ordering():
    """3 x 3 = 3bar + 6 forced. h(3bar) = 2/9 < h(6) = 5/9, so the
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
    """Under the normal-ordering benchmark, m1 = 0 plus measured
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


def test_bell_fork():
    """The local realist reading (outcomes as assigned +/-1 facts) is
    capped at 2. The framework's local non-realist reading (the
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
    constants are stored nowhere. They are where the recursion lands.
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




def test_alpha_pins_higgs_coordinate_and_G(res):
    """The two alpha_G2(M_Pl) pins stated in root.py, measured.
    (1) the transported Higgs coordinate responds to the bridge boundary.
    (2) G_ind: xi_bridge scales as alpha; +10% alpha shifts the
    heat-kernel sum Sigma by 0.64%, against G_ind/G_N closed at 8e-8."""
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

    # (2) the Sigma shift from xi_bridge under +10% alpha
    xi = float(rt.alpha_G2_Pl) * 0.5 * 1.0
    d_sigma = g["N_bridge"] * (-0.1 * xi) / g["Sigma_total_UV"]
    assert 0.005 < abs(d_sigma) < 0.008   # measured 0.64%
    assert abs(g["G_ratio_canonical"] - 1.0) < 1e-7


# ═══════════════════════════════════════════════════════════════════════
#  The derivation chain, asserted against the engine
#
#  These are the identities the chain stands on: the dimensional
#  transmutation actions and their two face readings, the depth-4
#  response weights, the per-mode condensate reading behind the
#  alpha_s threshold, and the CKM lambda form.  Each is computed in
#  the module that owns it; each is read back here so a refactor
#  cannot quietly move one.
# ═══════════════════════════════════════════════════════════════════════

def test_transmutation_actions_and_their_two_faces():
    """S_lepton = 2pi/(b0 alpha) on the F4 face, transported to the G2
    face by the h-dual metric: b_face = b0 h(F4)/h(G2) = 24, which is
    also the embedding invariant h(E8) - C2(26).  Hence S_baryo =
    (4/9) S_lepton = 2 pi^2, read at web depth 4."""
    import root

    assert root.b0_G2 == Fraction(32, 3)
    assert root.S_lepton == pytest.approx(
        2.0 * math.pi / (float(root.b0_G2) * root.alpha_G2_Pl), abs=1e-12)
    assert root.S_lepton == pytest.approx(
        9.0 * math.pi**2 / 2.0, abs=1e-12)
    assert root.b_baryo_face == root.hv_E8 - root.C2_26 == 24
    assert root.S_baryo == pytest.approx(
        (root.hv_G2 / root.hv_F4) * root.S_lepton, abs=1e-12)
    assert root.S_baryo == pytest.approx(2 * math.pi**2, abs=1e-12)
    assert root.WEB["S_baryo"].web_depth == 4
    assert root.alpha_G2_Pl == pytest.approx(1.0 / (24.0 * math.pi),
                                             abs=1e-15)
    assert root.delta_S == pytest.approx(15.0 / 512.0, abs=1e-15)


def test_depth4_response_weights():
    """The depth-3 e<->q carrier read through each node's own channel:
    action -> charge_trace * C2(26) = 16, amplitude -> W(fund) * dim
    J3(O) = 27/2.  The Albert lift is unital, so the junction
    normalization is one and no third weight enters."""
    import root

    assert float(root.charge_trace * root.C2_26) == 16.0
    assert root.d11**3 / root.d10 == 13.5


def test_condensate_map_is_the_per_mode_reading(res):
    """The G2 breaking scale is the connector's PER-MODE reading of the
    electroweak condensate: v7^2 = E[v^2] v_EW^2 with E[v^2] = 1/2, so
    M_V = g v_EW sqrt(C2(7)/6) with g = pi/(2 sqrt(2))."""
    import root

    v_ew = root.M_Pl_GeV * math.exp(-root.WEB.state["S_quark"])
    m_v = (math.pi / (2.0 * math.sqrt(2.0))) * v_ew * math.sqrt(
        float(root.E_v2) * float(root.C2_7) / 6.0)
    assert res["c"]["M_V_GeV"] == pytest.approx(m_v, abs=1e-2)
    assert float(root.E_v2) == 0.5


def test_ckm_lambda_is_the_base_plus_echo_form():
    """lambda = h10 (1 + h11/K^2) = (2/9)(73/72) = 73/324 exactly: the
    fundamental Sugawara weight carrying the SAME quadratic-altitude
    echo root.py applies to bridge^2."""
    import root

    assert 1 + root.h11 / root.K**2 == Fraction(73, 72)
    assert root.h10 * (1 + root.h11 / root.K**2) == Fraction(73, 324)
    assert abs(float(Fraction(73, 324)) - 0.2253086) < 1e-6


def test_lambda_echo_menu_unique_hit(res):
    """The committed backreaction h10(1 + h11/K^2) is the unique
    1 sigma hit over the engine's echo menu
    {h10, h11, Q0, 1} x {K, K^2, K^3}: the enumeration is executable,
    and the count is frozen here."""
    menu = res["mx"]["lambda_menu"]
    assert menu["n_candidates"] == 12
    assert menu["n_within_1sigma"] == 1
    assert menu["committed"] == ("h11", 2)
    assert menu["committed_is_unique_hit"]
    assert menu["hit_keys"] == [("h11", 2)]
    assert abs(menu["committed_pull"]) < 0.5
    assert menu["next_best_abs_pull"] > 1.9


def test_modular_phase_candidate_for_eta_bar():
    """The computed modular phase is pi/9; assigning it to eta_bar is IH."""
    import mixing
    import root

    wzw = mixing._compute_wzw()
    fund = wzw["IDX"][(1, 0)]
    phase = math.atan2(wzw["S_mat"][fund, fund].imag,
                       wzw["S_mat"][fund, fund].real)
    assert phase == pytest.approx(math.pi / root.d11**2, abs=1e-12)


def test_two_F4_singlets_available_for_RH_identification():
    """The 26 decomposition contains two singlets; identifying them as
    RH neutrinos gives a rank-at-most-two matrix only in minimal type-I."""
    import root

    assert root.n26 - 3 * (root.d11**2 - 1) == 2


def test_neutron_proton_splitting_qcd_part(res):
    """m_d - m_u = 5 m_e exactly, a datum the dictionary never
    consulted.  Lattice (BMW 2015): 2.52 +- 0.23 MeV."""
    m = res["m"]
    assert m["m_d"] - m["m_u"] == pytest.approx(5 * m["m_e"], abs=5e-3)


def test_confinement_scale_is_the_constituent_scale():
    """Lambda_conf = C_EFF M_Pl exp(-S_lepton) = 314.2 MeV, the
    constituent-quark scale to ~1%."""
    import root

    lam_conf = root.C_EFF * root.M_Pl_GeV * math.exp(-root.S_lepton)
    assert lam_conf * 1e3 == pytest.approx(314.2, abs=0.1)


def test_no_higgs_mass_leaves_the_chain(res):
    """The mixed-scheme 125.2965 reading is not canonical: no Higgs
    mass prediction leaves the derivation."""
    assert "mH_pred" not in res["g"]
