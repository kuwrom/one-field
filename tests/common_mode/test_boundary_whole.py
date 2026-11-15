"""Freezes tests/probes/boundary_whole.py: the boundary reading c as the common mode reading no share."""

import math

import pytest


def _silent(*args, **kwargs):
    pass


def test_the_common_mode_reads_each_sum_whole():
    """The pointed probability 1/12, the local response 1/(2 pi) and
    the per channel moment 1/2 are shares of sums that are each one,
    and the engine's alpha_G2, xi and dark ratio are products of
    shares. Of the twelve readings of c built from them, four pass
    gates S and M with no measured G, and only c = 1 passes S under all
    six conventions.  c = 2/B^2 moves with the harmonic, and the lepton
    ratio holds it to 9e-10.  With c = 1 read as the common mode's depth
    four count, the anchor's G is 54.6 ppm below the recommended
    value."""
    from probes import boundary_whole

    r = boundary_whole.run(report=_silent)
    assert r["alpha_G2"] == pytest.approx(1 / (24 * math.pi), abs=1e-15)
    assert r["dark_ratio"] == pytest.approx(2 * math.pi - 1, abs=1e-12)
    assert r["coefficient"] == pytest.approx(1.0, abs=1e-12) and r["c"] == pytest.approx(1.0, abs=1e-12)
    assert len(r["menu"]) == 12
    assert r["passing"] == pytest.approx([1.0, 1.1893, 2.0, 2.3786], abs=1e-3)
    assert [k for k, n in r["conventions"].items() if n == 6] == ["whole/whole/whole"]
    assert r["dark_readings"]["local alone"] == 0.0 and r["dark_readings"]["remainder alone"] == math.inf
    assert r["compression_error"] < 1e-12
    assert r["c_moves"] == pytest.approx(-1e-3, abs=1e-5) and r["ratio_slope"] == pytest.approx(24.0, abs=1.0)
    assert r["c_held"] < 2e-9 and r["Q"] == pytest.approx(2 / 3)
    assert r["common_mode_coefficient"] == pytest.approx(1.0) and abs(r["mean_over_frame"]) < 2e-3
    assert r["g_ppm"] == pytest.approx(-54.6, abs=0.1) and r["g_sigma"] == pytest.approx(-2.43, abs=0.01)


def test_the_channel_map_is_the_moments_whole():
    """On the searchlight's circulant the common mode's count is
    (B/A)^2/2 = 3 Q - 1 for every B/A, so the channel map sqrt d10, the
    Koide value 2/3 and the whole second moment are one statement.
    Without data, the amplitude gate's window puts the count between
    0.852 and 1.086.  G holds B/A to 1.1e-5 and m_mu/m_e to 5e-10."""
    from probes import boundary_whole

    r = boundary_whole.run(report=_silent)
    assert all(v == pytest.approx(ba * ba / 2, abs=1e-12) for ba, v in r["channel_map"].items())
    assert r["ba_held"]["amplitude gate"] == pytest.approx(0.059, abs=0.001)
    assert r["ba_held"]["G"] == pytest.approx(1.12e-5, abs=1e-7) and r["ba_held"]["m_mu/m_e"] < 1e-9
    assert r["classical_count"] == pytest.approx(0.25, abs=1e-12)
