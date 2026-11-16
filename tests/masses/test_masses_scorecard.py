"""Scorecard tests of the masses branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Frozen mass table  (MeV)
# ═══════════════════════════════════════════════════════════════════════


FROZEN_MASSES_MEV = {
    'm_e': 0.51099895069,        # anchor, exact by construction
    'm_mu': 105.6583758,
    'm_tau': 1776.9092814,
    'm_u': 2.157551125,
    'm_d': 4.712545879,
    'm_s': 93.8395306,           # depth-2 bridge value
    'm_c': 1273.770419,
    'm_b': 4193.80675,
    'm_t': 172508.2761,
}


# ═══════════════════════════════════════════════════════════════════════
#  Quark masses  (6 predictions)
# ═══════════════════════════════════════════════════════════════════════

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


def test_strange_koide_is_649_972(res):
    """Q(s,c,b) = 2/3 + h(1,0)/K^3 = 649/972 (Rivero fundamental channel).

    masses.py solves the strange root from the ledger value.
    """
    import root
    K = 6
    Q = Fraction(2, 3) + Fraction(2, 9) / Fraction(K ** 3, 1)
    assert Q == Fraction(649, 972)
    assert root.WEB["Q_scb"].value() == pytest.approx(649 / 972, rel=1e-12)


def test_albert_bridge_factor_is_32_over_27():
    """bridge^2 = Q_0^2 * d_10^3 / d_11 = (4/9)(8/3) = 32/27 is the ledger
    base. The shipped value carries the 73/72 echo: 292/243."""
    import root
    assert float(root.Q0**2 * root.charge_trace) == 32 / 27
    assert root.WEB.state["bridge_sq"] == pytest.approx(292 / 243, rel=1e-12)


def test_lepton_koide_is_two_thirds(res):
    """Q_lep = 1/3 + (B/A)^2 / 6 = 2/3 from |B/A|^2 = 2."""
    Q = Fraction(1, 3) + Fraction(2, 6)
    assert Q == Fraction(2, 3)
    m = res["m"]
    me, mmu, mtau = m["m_e"], m["m_mu"], m["m_tau"]
    Q_num = (me + mmu + mtau) / (math.sqrt(me) + math.sqrt(mmu) + math.sqrt(mtau)) ** 2
    assert math.isclose(Q_num, 2 / 3, rel_tol=1e-5)


# ═══════════════════════════════════════════════════════════════════════
#  Generation word lemma
# ═══════════════════════════════════════════════════════════════════════

def test_generation_words_are_walk_counts():
    """base(n) = (4, 12, 97) as boundary walk counts. Z2 split u/d."""
    import words
    w = words.derive()
    assert w["base"] == (4, 12, 97)
    assert w["base_down_seed"] == 9          # d = aa (Z2-swapped seed)
    assert w["FP_C0"] == 12                   # neutral lane = base(2)
    # mixed Z2-symmetric terminal words are excluded (data: 64/139 GeV)
    assert set(w["excluded"].values()) == {36, 78}


# ═══════════════════════════════════════════════════════════════════════
#  Vertex composition rule + anchor inversion (m_e → M_Pl → G)
# ═══════════════════════════════════════════════════════════════════════

def test_electron_anchor_exact(res):
    """m_e reproduces the sole dimensional input exactly."""
    from root import M_E_ANCHOR_MEV
    assert abs(res["m"]["m_e"] / M_E_ANCHOR_MEV - 1.0) < 1e-12


def test_muon_mass_predicted_ppb(res):
    """m_mu is PREDICTED below the ppm level, within 1 sigma of CODATA."""
    m_mu_codata, sigma = 105.6583755, 2.3e-6
    assert abs(res["m"]["m_mu"] - m_mu_codata) / sigma < 1.0


def test_tau_mass_prediction_band(res):
    """The m_tau prediction lies inside the PDG band."""
    assert abs(res["m"]["m_tau"] - 1776.93) / 0.09 < 1.0


# ═══════════════════════════════════════════════════════════════════════
#  The canonical freeze: the full value table, pinned.
#  An edit that moves a canonical number past these tolerances is
#  either a derivation change (update it deliberately, with provenance)
#  or a regression (fix it).
# ═══════════════════════════════════════════════════════════════════════

def test_mass_coordinate_invariants(res):
    """Scheme invariant light quark ratios. These predictions need no
    choice of mass coordinate, because QCD running cancels."""
    ms = res["m"]
    assert abs(ms["mu_over_md"] - 38.0 / 83.0) < 1e-12   # exact
    assert abs(ms["ms_over_mud"] - 27.3183) < 1e-3
    assert abs(ms["Q_ellipse"] - 22.383) < 1e-2


@pytest.mark.parametrize("key,val", sorted(FROZEN_MASSES_MEV.items()))
def test_frozen_mass_table(res, key, val):
    """Every mass the structure produces, frozen to 10 digits."""
    assert abs(res["m"][key] / val - 1.0) < 1e-9


def test_np_splitting_qcd_part(res):
    """m_d - m_u = 5 m_e exactly. This is compared with the lattice QCD
    part of m_n - m_p (BMW, Science 347:1452 (2015): 2.52 +/- 0.23 MeV)."""
    m = res["m"]
    dmq = m["m_d"] - m["m_u"]
    assert abs(dmq / (5.0 * m["m_e"]) - 1.0) < 1e-12
    pull = (dmq - 2.52) / 0.23
    assert abs(pull) < 2.0          # +0.15 sigma


def test_diquark_channel_ordering():
    """3 x 3 = 3bar + 6 is forced. Since h(3bar) = 2/9 < h(6) = 5/9,
    the nucleon (good diquark) is forced to sit below the Delta."""
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


def test_depth4_response_weights():
    """The depth-3 e<->q carrier, read through each node's own channel,
    gives two weights: action -> charge_trace * C2(26) = 16 and
    amplitude -> W(fund) * dim J3(O) = 27/2. The Albert lift is unital,
    so the junction normalization is one and no third weight enters."""
    import root

    assert float(root.charge_trace * root.C2_26) == 16.0
    assert root.d11**3 / root.d10 == 13.5


def test_confinement_scale_is_the_constituent_scale():
    """Lambda_conf = C_EFF M_Pl exp(-S_lepton) = 314.2 MeV, the
    constituent quark scale to ~1%."""
    import root

    lam_conf = root.C_EFF * root.M_Pl_GeV * math.exp(-root.S_lepton)
    assert lam_conf * 1e3 == pytest.approx(314.2, abs=0.1)
