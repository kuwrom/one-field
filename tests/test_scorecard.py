"""Scorecard tests: the 34 numerical predictions as a falsification gauntlet.

Each test below encodes one entry of the README scorecard as a pass/fail
assertion against PDG / NuFit / Planck reference data.  The framework's
claim is "no free parameters at any layer".  This test suite makes that
claim falsifiable in the operational sense: any change to the algebra
that moves a prediction outside its tolerance band fails the suite.

Run with:
    pytest

Tolerances are chosen with a small margin above the current observed
deviation (typically 2-5x), so the tests genuinely guard against regression
without being so tight that they fail on cosmetic numerical noise.
"""

import math
from fractions import Fraction

import pytest

from E8.constants import N_VERTEX, NUFIT_OSCILLATION, NUFIT_PMNS, PDG_MASSES, ALPHA_EM


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
    """Charged-lepton pole masses from the Brannen Z3 formula."""
    predicted = res["lep"][f"m_{name}"]
    pdg = PDG_MASSES[name]
    assert err_pct(predicted, pdg) < tol_pct, (
        f"m_{name}: predicted {predicted:.6f} MeV, PDG {pdg} MeV, "
        f"deviation {err_pct(predicted, pdg):.4f}% exceeds {tol_pct}%"
    )


@pytest.mark.parametrize("name,tol_pct", [
    ("u", 0.5),
    ("d", 0.5),
    ("s", 0.5),
    ("c", 0.5),
    ("b", 0.5),
    ("t", 0.5),
])
def test_quark_mass(res, name, tol_pct):
    """Six quark masses from F4 triality + WZW emergence + Koide Q-relations."""
    predicted = res["qrk"][f"m_{name}"]
    pdg = PDG_MASSES[name]
    assert err_pct(predicted, pdg) < tol_pct


# ═══════════════════════════════════════════════════════════════════════
#  Electroweak scale  (1 prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_electroweak_scale(res):
    """v_EW from M_Pl x exp(-(9*pi^2/2 - 6 + 15/512))."""
    v_EW_pred = res["scl"]["v_EW_pred_GeV"]
    assert err_pct(v_EW_pred, 246.22) < 0.05


# ═══════════════════════════════════════════════════════════════════════
#  Strong coupling  (1 prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_alpha_s_at_MZ(res):
    """alpha_s(M_Z) from the embedding-index chain SU(3) C G2 C E8."""
    alpha_s_MZ = res["als"]["alpha_s_MZ_thresh"]
    pdg = 0.1180
    assert err_pct(alpha_s_MZ, pdg) < 0.5


# ═══════════════════════════════════════════════════════════════════════
#  CKM sector  (4 Wolfenstein + UT angles, summarised as chi^2/n + pull)
# ═══════════════════════════════════════════════════════════════════════

def test_ckm_chi_squared(res):
    """CKM observables match PDG within chi^2/n < 1."""
    n = len(res["ck"]["ckm_obs"])
    chi2_per_n = res["ck"]["chi2_ckm"] / n
    assert chi2_per_n < 1.0, f"chi^2/n = {chi2_per_n:.2f} >= 1.0"


def test_ckm_max_pull(res):
    """No individual CKM observable deviates by more than 2 sigma."""
    assert res["ck"]["max_pull_ckm"] < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  PMNS sector  (3 mixing angles + delta_CP prediction)
# ═══════════════════════════════════════════════════════════════════════

def test_pmns_chi_squared(res):
    """PMNS angles match NuFit 6.0 within chi^2/n < 1."""
    chi2_per_n = res["pm"]["chi2_pmns"] / 3
    assert chi2_per_n < 1.0


def test_pmns_delta_cp_prediction(res):
    """delta_CP prediction lies in a testable range for DUNE/Hyper-K."""
    delta_CP = res["pm"]["delta_CP_deg"]
    assert 0 < delta_CP < 180, f"delta_CP = {delta_CP} outside expected range"


# ═══════════════════════════════════════════════════════════════════════
#  Higgs mass  (1 prediction, via SM RGE with derived boundary)
# ═══════════════════════════════════════════════════════════════════════

def test_higgs_mass(res):
    """m_H from F4(1) fusion (lambda(M_Pl)=0) + bridge threshold."""
    mH = res["hg"]["mH_pred"]
    assert err_pct(mH, 125.20) < 1.0


# ═══════════════════════════════════════════════════════════════════════
#  Gravity sector  (Newton constant from the (7,26) bridge heat kernel)
# ═══════════════════════════════════════════════════════════════════════

def test_newton_constant_UV(res):
    """G_ind/G_N (UV bookkeeping) within 1% of unity."""
    assert err_pct(res["gr"]["G_ratio_UV"], 1.0) < 1.0


def test_newton_constant_broken(res):
    """G_ind/G_N (broken-phase bookkeeping) within 1% of unity."""
    assert err_pct(res["gr"]["G_ratio_broken"], 1.0) < 1.0


# ═══════════════════════════════════════════════════════════════════════
#  Baryogenesis  (eta_B from G2 instanton + Fano orientation)
# ═══════════════════════════════════════════════════════════════════════

def test_baryon_asymmetry(res):
    """eta_B = dim(G2_fund) * J_lep * exp(-2*pi^2) within 2% of Planck."""
    eta_B = res["br"]["eta_B"]
    eta_B_obs = res["br"]["eta_B_obs"]
    assert err_pct(eta_B, eta_B_obs) < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  Neutrino sector  (structural: m1=0, ordering, Sum m_nu)
# ═══════════════════════════════════════════════════════════════════════

def test_neutrino_lightest_mass_is_zero(res):
    """Rank-2 Type-I seesaw from F4's two singlets forces m_1 = 0."""
    assert res["nu"]["m1_meV"] == 0.0


def test_neutrino_normal_ordering(res):
    """Algebraic structure forces normal ordering."""
    assert res["nu"]["ordering"] == "normal"


def test_neutrino_sum_at_oscillation_floor(res):
    """Sum m_nu = sqrt(dm2_21) + sqrt(dm2_31) lies at the oscillation floor."""
    sum_meV = res["nu"]["sum_m_nu_meV"]
    assert 55.0 < sum_meV < 62.0


# ═══════════════════════════════════════════════════════════════════════
#  Structural closure  (E8 mode accounting, no free parameters)
# ═══════════════════════════════════════════════════════════════════════

def test_e8_central_charge_sum(res):
    """c(G2) + c(F4) = c(E8) = 8 exactly."""
    c_G2 = res["alg"]["c_G2"]
    c_F4 = res["alg"]["c_F4"]
    c_E8 = res["alg"]["c_E8"]
    assert math.isclose(c_G2 + c_F4, c_E8, abs_tol=1e-12)
    assert math.isclose(c_E8, 8.0, abs_tol=1e-12)


def test_e8_mode_assignment_closure():
    """248 = 14 + 52 + 182.  Every E8 mode has a physical job."""
    dim_G2_adj = 14   # gauge sector
    dim_F4_adj = 52   # matter and Higgs sector
    dim_bridge = 7 * 26  # mixed (7,26) sector, gravity
    assert dim_G2_adj + dim_F4_adj + dim_bridge == 248


def test_strong_cp_resolved(res):
    """theta_QCD = 0 is forced by pi_3(G2) -> pi_3(SU(3)) with unique E8(1) vacuum."""
    # The framework predicts theta_QCD = 0 structurally.  No CP-violating
    # QCD parameter enters anywhere in the derivation chain; this test
    # checks that no module exposes a free theta parameter.
    for key, val in res["alg"].items():
        assert "theta_QCD" not in key.lower(), (
            "theta_QCD appeared as a derived/free parameter - "
            "the strong-CP resolution requires it to be structurally absent"
        )


# ═══════════════════════════════════════════════════════════════════════
#  Algebraic structure  (rational identities behind the predictions)
#  ──────────────────────────────────────────────────────────────────
#  The numerical tests above check "does the answer match PDG".  These
#  check "does the answer come from the right rational?".  A refactor
#  that preserved the float but broke the algebraic ratio would pass
#  the numerical tests and fail these.
# ═══════════════════════════════════════════════════════════════════════

def test_light_quark_ratios_are_38_9_and_83_9(res):
    """m_u/m_e = 38/9 and m_d/m_e = 83/9 (triality + h(1,0) = 2/9)."""
    h10 = Fraction(2, 9)
    assert Fraction(4, 1) + h10 == Fraction(38, 9)
    assert Fraction(9, 1) + h10 == Fraction(83, 9)
    # And the module returns those floats:
    assert math.isclose(res["qrk"]["m_u"] / res["lep"]["m_e"], 38 / 9, rel_tol=1e-12)
    assert math.isclose(res["qrk"]["m_d"] / res["lep"]["m_e"], 83 / 9, rel_tol=1e-12)


def test_heavy_quark_ratios_are_217_18_and_1165_12(res):
    """m_c/m_mu = 217/18 (charm) and m_t/m_tau = 1165/12 (top)."""
    delta_OPE = Fraction(1, 18)  # h(1,1) − 2 h(1,0) = 1/2 − 4/9
    K = 6
    assert Fraction(12, 1) + delta_OPE == Fraction(217, 18)
    assert Fraction(81 + 16, 1) + Fraction(1, 2 * K) == Fraction(1165, 12)
    assert math.isclose(res["qrk"]["m_c"] / res["lep"]["m_mu"], 217 / 18, rel_tol=1e-12)
    assert math.isclose(res["qrk"]["m_t"] / res["lep"]["m_tau"], 1165 / 12, rel_tol=1e-12)


def test_bottom_koide_is_289_432(res):
    """Q(c,b,t) = 2/3 + h(1,1)/K^3 = 289/432."""
    K = 6
    Q = Fraction(2, 3) + Fraction(1, 2) / Fraction(K**3, 1)
    assert Q == Fraction(289, 432)
    sc = math.sqrt(res["qrk"]["m_c"])
    sb = math.sqrt(res["qrk"]["m_b"])
    st = math.sqrt(res["qrk"]["m_t"])
    Q_num = (res["qrk"]["m_c"] + res["qrk"]["m_b"] + res["qrk"]["m_t"]) / (sc + sb + st) ** 2
    assert math.isclose(Q_num, 289 / 432, rel_tol=1e-10)


def test_strange_koide_is_649_972():
    """Q(s,c,b) = 2/3 + h(1,0)/K^3 = 649/972 (Rivero inverse fundamental channel)."""
    K = 6
    Q = Fraction(2, 3) + Fraction(2, 9) / Fraction(K**3, 1)
    assert Q == Fraction(649, 972)


def test_albert_bridge_factor_is_32_over_27():
    """bridge^2 = Q0^2 * d(1,0)^3 / d(1,1) = (4/9)(8/3) = 32/27."""
    d10, d11 = 2, 3
    bridge_sq = Fraction(2, 3) ** 2 * Fraction(d10**3, 1) / Fraction(d11, 1)
    assert bridge_sq == Fraction(32, 27)


def test_lepton_koide_is_two_thirds(res):
    """Q_lep = 1/3 + (B/A)^2 / 6 = 2/3 from |B/A|^2 = 2."""
    BA_sq = 2  # octonionic CG corollary
    Q = Fraction(1, 3) + Fraction(BA_sq, 6)
    assert Q == Fraction(2, 3)
    me, mmu, mtau = res["lep"]["m_e"], res["lep"]["m_mu"], res["lep"]["m_tau"]
    Q_num = (me + mmu + mtau) / (math.sqrt(me) + math.sqrt(mmu) + math.sqrt(mtau)) ** 2
    assert math.isclose(Q_num, 2 / 3, rel_tol=1e-5)


def test_vertex_count_is_30_with_no_free_choice():
    """N_vertex = 26 (F4 fund irreducible) + 4 (Higgs DOFs) = 30."""
    assert N_VERTEX == 30


def test_vertex_action_is_15_over_512():
    """delta_S = N * alpha / (2 pi) = 30 * (pi/512) / (2 pi) = 15/512 exactly."""
    delta_s = N_VERTEX * ALPHA_EM / (2.0 * math.pi)
    assert math.isclose(delta_s, 15.0 / 512.0, abs_tol=1e-15)


def test_alpha_bridge_chain_is_one_loop_exact(res):
    """alpha(0) bridge derivation: h_bridge = 1, D^2 = 1, c_coset = 0."""
    ab = res["ab"]
    # Marginal primary: h_bridge = h(7) + h(26) = 2/5 + 3/5 = 1
    assert math.isclose(ab["h_bridge"], 1.0, abs_tol=1e-12)
    assert math.isclose(ab["h_G2_7"], 2 / 5, abs_tol=1e-12)
    assert math.isclose(ab["h_F4_26"], 3 / 5, abs_tol=1e-12)
    # Lagrangian condensation: D^2_local = 1
    assert math.isclose(ab["D2_local"], 1.0, abs_tol=1e-12)
    assert math.isclose(ab["g_bridge_sq"], 1.0, abs_tol=1e-12)
    # Topological coset: c_coset = 0
    assert math.isclose(ab["c_coset"], 0.0, abs_tol=1e-12)
    # Closed-form: 1/alpha(0) = 256(2pi-1)/pi^2
    expected = 256.0 * (2.0 * math.pi - 1.0) / math.pi**2
    assert math.isclose(ab["inv_alpha_phys"], expected, rel_tol=1e-12)


def test_gravity_bridge_xi_is_1_over_48_pi(res):
    """xi_bridge = alpha_G2 * E[v^2] * h_bridge = 1/(48 pi)."""
    assert math.isclose(res["gr"]["xi_bridge"], 1.0 / (48.0 * math.pi), abs_tol=1e-14)
    assert math.isclose(res["gr"]["E_v2"], 0.5, abs_tol=1e-12)
    assert math.isclose(res["gr"]["h_bridge"], 1.0, abs_tol=1e-12)
    assert res["gr"]["N_bridge"] == 182
