"""Freezes tests/probes/g_measurements.py: the anchor's G against the
sixteen measurements behind the recommended value, and the precision
block's joint statistic."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_g_measurements_anchor_sits_in_the_scatter_below_its_precise_half(res):
    """The sixteen inputs of Table XXX (arXiv:2409.03787) give the mean
    6.674300(38) with chi-square 195.7 on 15, and 3.9 times that
    uncertainty is the recommended 6.674 30(15).  The anchor's 6.673936
    lies 9.6 unexpanded and 2.47 expanded standard deviations below the
    recomputed mean, and 2.43 below the recommended value.  Seven values
    lie below it and nine above.  It lies below each of the five most
    precise inputs by 2.6 to 7.1 of their own standard deviations.  The
    block G, m_mu, m_tau, v_EW has chi-square 6.0 on 4, of which G
    carries 0.99.  Under an added spread, the unweighted mean, the
    median and each of the two published consensus values, the anchor
    lies within one standard deviation.  The 2.43 belongs to the
    multiplicative rule.  The anchor read through web depth one and the
    boundary read as 1/(4 pi G) fail the comparison."""
    from probes import g_measurements as gm

    out = gm.run(report=_silent, build=False)
    assert (out["mean"], out["mean_unc"]) == pytest.approx((6.674300, 3.8e-5), abs=2e-6)
    assert out["chi2_set"] == pytest.approx(195.7, abs=0.1) and out["dof"] == 15
    assert out["raw"] == pytest.approx(-9.63, abs=0.02) and out["expanded"] == pytest.approx(-2.47, abs=0.01)
    assert (out["below"], out["above"]) == (7, 9) and out["neighbours"] == ("MSL-03", "LANL-97")
    assert out["pulls"]["JILA-18"] == pytest.approx(-5.34, abs=0.01)
    assert out["pulls"]["BIPM-14"] == pytest.approx(10.03, abs=0.01)
    assert sorted(out["precise"]) == ["HUSTA-18", "HUSTT-18", "UCI-14", "UWash-00", "UZur-06"]
    assert min(out["precise"].values()) == pytest.approx(2.62, abs=0.01)
    assert max(out["precise"].values()) == pytest.approx(7.12, abs=0.01)
    b = out["block"]
    assert [b[k] for k in ("G", "m_mu", "m_tau", "v_EW")] == pytest.approx([-2.43, 0.12, -0.23, 0.08], abs=0.01)
    assert out["chi2_block"] == pytest.approx(5.97, abs=0.02)
    assert out["p_block"] == pytest.approx(0.20, abs=0.01) and out["p_worst"] == pytest.approx(0.06, abs=0.01)
    rules = out["rules"]
    assert rules["added spread, chi-square = n - 1"][:2] == pytest.approx((6.67389, 0.00026), abs=1e-5)
    assert rules["added spread, moment estimate"][:2] == pytest.approx((6.67399, 0.00016), abs=1e-5)
    assert rules["median"][:2] == pytest.approx((6.67408, 0.00023), abs=1e-5)
    assert [round(v[2], 1) for v in rules.values()] == [0.2, -0.4, 0.6, -0.6]
    assert out["published"]["Merkatas et al."] == pytest.approx(-0.60, abs=0.01)
    assert abs(out["published"]["Rinaldi et al."]) < 0.05
    assert out["published"]["Muravyov et al."] == pytest.approx(-11.46, abs=0.01)
    assert out["alpha"]["rubidium recoil"] == pytest.approx(-14.2, abs=0.1)
    assert out["alpha"]["caesium recoil"] == pytest.approx(0.13, abs=0.01)


def test_g_measurements_comparison_rejects_other_readings():
    """A value passes when it lies within three expanded standard
    deviations of the mean and inside the set's scatter.  The recommended
    value passes; a G 2280 ppm high (the anchor at web depth one) and a
    G smaller by 4 pi do not."""
    from probes import g_measurements as gm

    assert gm.compare(6.673936)["passes"] and gm.compare(gm.RECOMMENDED[0])["passes"]
    assert not gm.compare(6.673936 * (1 + 2280e-6))["passes"]
    assert not gm.compare(6.673936 / 12.566)["passes"]
    assert gm.compare(6.6770)["below"] == 16
