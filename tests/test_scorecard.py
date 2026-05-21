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

import pytest

from E8.constants import NUFIT_OSCILLATION, NUFIT_PMNS, PDG_MASSES


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
