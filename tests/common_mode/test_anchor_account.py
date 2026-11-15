"""Freezes tests/audit/anchor_account.py: the anchor relation held as one
object, its readers computed, and the homes of G's residual."""

import math

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def books():
    import anchor_account
    return anchor_account, anchor_account.run(report=_silent)


def test_anchor_account_nothing_off_the_books(books):
    """Ten entries reproduce the engine's M_Pl, m_mu, m_tau and v_EW to
    1e-12, and the hand written inversions agree with them."""
    _, out = books
    assert out["entries"] == 10 and out["completeness"] < 1e-12 and out["copies"] < 1e-12
    assert out["pulls"]["G"] == pytest.approx(-2.43, abs=0.01)


def test_anchor_account_residual_has_one_home_read_by_g_alone(books):
    """Each entry's readers are found by perturbing it.  c moves G
    alone.  S0 moves G and 1/alpha(0), whose base is the endpoint of the
    lever arm that S0 measures.  Only for c does a shift that zeroes G's
    pull leave the other observables untouched.  1/alpha(0) excludes S0
    at twenty three thousand standard deviations.  G_F excludes C_EFF and
    c2 at 106.  The lepton ratios exclude theta, B/A and w4 at a
    thousand or more."""
    _, out = books
    assert out["readers"]["S0"] == ["G", "1/alpha"] and out["readers"]["c"] == ["G"]
    assert out["read_by_G_alone"] == ["c"] and out["homes"] == ["c"]
    assert out["cost"]["S0"][1] == pytest.approx(23088, rel=0.01)
    assert out["cost"]["c"][0] == pytest.approx(-5.46e-5, rel=0.01) and out["cost"]["c"][1] < 0.3
    assert out["readers"]["w16"] == ["v_EW"] and "m_mu" in out["readers"]["w4"]
    assert out["cost"]["S0"][0] == pytest.approx(-2.73e-5, rel=0.01)
    assert out["pulls"]["1/alpha"] == pytest.approx(0.13, abs=0.01)
    assert out["cost"]["C_EFF"][1] == pytest.approx(106, abs=1) and out["cost"]["c2"][1] == pytest.approx(106, abs=1)
    assert out["cost"]["theta"][1] > 1000 and out["cost"]["w4"][1] == pytest.approx(2127, rel=0.01)
    assert math.isinf(out["tolerance"]["c"]) and out["tolerance"]["S0"] < 1e-8 and out["tolerance"]["theta"] < 1e-9
    assert out["w16_band"] == pytest.approx(0.38, abs=0.01)


def test_anchor_account_refuses_what_is_not_booked(books):
    """A shift of a name that is not an entry is refused, and the
    comparison constant's rounding is stated: 0.105 ppm, the engine's G
    from primaries 6.6739371e-11 at -2.42."""
    module, out = books
    with pytest.raises(AssertionError):
        module.Account().observe(w_frame=1.0)
    assert out["rounding_ppm"] == pytest.approx(-0.105, abs=0.002)
    assert out["G_direct"] == pytest.approx(6.6739371e-11, abs=1e-17)
    assert out["pull_direct"] == pytest.approx(-2.419, abs=0.002)
