"""Freezes tests/audit/carrier_ledger.py: every node of the web against the closed carrier."""

import pytest

import carrier_ledger


@pytest.fixture(scope="module")
def ledger():
    return carrier_ledger.run(report=lambda *a: None)


def test_six_read_ten_follow_eleven_read_none(ledger):
    """The readers are the coupling, the quark action, the three lepton
    rates and the common mode.  The followers reach a reader through
    the web's own reads.  The anchor holds m_e fixed, and with it the
    two quarks read from m_e."""
    assert set(ledger["readers"]) == {"inv_alpha", "S_quark", "m_e", "m_mu", "m_tau", "common_mode"}
    assert ledger["followers"] == ["M_Pl_MeV", "lepton_EM", "m_b", "m_c", "m_d", "m_s", "m_t", "m_u",
                                   "sin2W", "v_EW_GeV"]
    assert ledger["still"] == ["G_ratio", "Q_cbt", "Q_scb", "S_G2", "S_lepton", "Sigma_vent(UV)",
                               "Sigma_vent(broken)", "a1_bridge", "bridge_sq", "dark_ratio", "lambda_MPl"]
    assert set(ledger["moved"]) == (set(ledger["readers"]) | set(ledger["followers"])) - {"m_e", "m_u", "m_d"}
    assert ledger["response"]["common_mode"] == pytest.approx(-7.502, abs=2e-3)
    # v_EW follows the action's 16 and the ruler's 0.545, the electron's own term
    assert ledger["response"]["v_EW_GeV"] == pytest.approx(16.0 + 0.5447, abs=0.01)
    assert ledger["response"]["M_Pl_MeV"] == pytest.approx(0.5447, abs=0.005)


def test_the_signs_and_the_weight_is_the_nodes_own(ledger):
    """On the ledgers of the coupling, the action and the rates, the term
    is negative in all five.  The common mode carries no weight: its two
    edges sum to one on the searchlight's circulant, and the carrier
    moves its count by -7.50 L through the three lepton terms.  A
    lepton's weight is half of the whole traversal 27, times its gap."""
    assert len(ledger["signs"]) == 5 and all(v < 0 for v in ledger["signs"].values())
    assert ledger["own"] == pytest.approx(1.0, abs=1e-14) and ledger["traversal"] == 27
    assert ledger["readers"]["common_mode"] == pytest.approx(-7.502, abs=1e-3)
    assert ledger["readers"]["S_quark"] == pytest.approx(-16.0, abs=1e-12)


def test_what_the_data_allow(ledger):
    """The allowed weights: 2.0000 +- 0.0001 by Cs recoil (1.9993 by Rb,
    1.9995 by the electron moment), 16.0 +- 0.2 by G_F, 13.50 +- 0.03
    by m_mu/m_e, and 33 +- 17 for the common mode's slope by G.  The
    count's -7.50 sits 2.43 standard deviations below the last.  A unit
    weight on the lepton action moves 1/alpha by over a thousand
    standard deviations.  On sin2W or the dark ratio it moves them by
    under a tenth of one."""
    a = ledger["allow"]
    assert a["1/alpha, Cs recoil"] == pytest.approx((2.0, 1.2e-4), abs=2e-5)
    assert a["1/alpha, Rb recoil"][0] == pytest.approx(1.9993, abs=1e-4)
    assert a["1/alpha, electron moment"][0] == pytest.approx(1.9995, abs=1e-4)
    assert a["S_quark, G_F"] == pytest.approx((15.99, 0.19), abs=0.01)
    assert a["leptons, m_mu/m_e"] == pytest.approx((13.504, 0.031), abs=2e-3)
    assert a["common mode, G"] == pytest.approx((33.0, 16.7), abs=0.1)
    slope = ledger["readers"]["common_mode"]
    assert (slope - a["common mode, G"][0]) / a["common mode, G"][1] == pytest.approx(-2.43, abs=0.01)
    u = ledger["unit"]
    assert u["S_lepton (1/alpha, Cs recoil)"] == pytest.approx(1140, rel=0.01)
    assert u["sin2W"] < 0.1 and u["dark_ratio"] < 0.1


def test_each_readers_sign_by_its_own_data(ledger):
    """The opposite sign misses Cs by 3.3e4 sigma, G_F by 168, m_mu/m_e by 884."""
    f = ledger["flip"]
    assert f["1/alpha, Cs recoil"] == pytest.approx(3.26e4, rel=0.01)
    assert f["S_quark, G_F"] == pytest.approx(-168, abs=2)
    assert f["leptons, m_mu/m_e"] == pytest.approx(-884, abs=5)
