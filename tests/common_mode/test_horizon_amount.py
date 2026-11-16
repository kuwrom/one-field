"""Freezes tests/probes/horizon_amount.py: gravity's sign and size read only the amount past the searchlight's horizon."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight():
    """The nameable rows give -4.5 at the covariant contact (-7.17 with
    the metric loop), which is repulsive.  Neither an octonion's worth of
    real scalars nor the named G2 singlets make them attractive.  A light
    amount of 80.7 to 392.5 real scalar units passes S and M at every
    contact with the loop (64.7 to 376.5 without), for any cutoff shape
    and placement.  The searchlight's 182 lies inside this window with a
    factor of about two to spare.  The label's coupling on the amount
    moves the edges by 4 percent.  With a mass, the amount keeps both
    gates below 0.672 M_b (0.574 M_b with the loop) under the least
    favourable monotone regulator.  Band pass shapes at the spectral
    amplitude need 0.05 to 0.10 M_b.  The labelled bridge keeps the
    anchor only above 0.999945 M_b, where the ledger's cutoff leaves it
    repulsive.  The one list's kernel is attractive and keeps the
    Einstein sign on the axis, with the limit 65/8.  A zero falls inside
    the boundary only below 139.3 units.  The label's coupling on the
    amount would put one at 1.295 M_b.  At the count the gates hold for
    -0.178 < xi < 0.107, and the zero stays outside the boundary for
    xi < 0.0151.  Of the five readings the list can give the bridge,
    only the amount keeps b0 and passes S and M."""
    from probes import boundary_content as bc
    from probes import horizon_amount as ha

    r = ha.run(report=_silent)
    assert ha.NAMED_SINGLETS == pytest.approx(bc.run(report=_silent)["stack"], abs=2e-3)
    assert r["named"] == pytest.approx({1.0: -4.5, 0.0: 13.5, -1.0: 31.5}, abs=1e-12)
    assert r["named_loop"] == pytest.approx(-7.1667, abs=1e-4) and r["octonion"] == pytest.approx(4 / 3)
    assert r["attract"] == pytest.approx(27.0, abs=1e-2) and r["attract_loop"] == pytest.approx(43.0, abs=1e-2)
    win = r["window"]
    assert win["no loop"] == pytest.approx((64.70, 376.49), abs=1e-2)
    assert win["loop"] == pytest.approx((80.70, 392.49), abs=1e-2)
    assert win["no loop, label's coupling"] == pytest.approx((67.38, 392.09), abs=1e-2)
    assert win["loop, label's coupling"] == pytest.approx((84.04, 408.75), abs=1e-2)
    assert r["margins"]["loop"] == pytest.approx((2.255, 2.157), abs=1e-3)
    assert r["margins"]["no loop"] == pytest.approx((2.813, 2.069), abs=1e-3)
    assert list(r["sigma"].values()) == pytest.approx([25.8333, 43.8333, 61.8333], abs=1e-4)
    assert r["t0"][1.0] == pytest.approx(4.1115, abs=1e-4)
    assert r["m_star"]["no loop"] == pytest.approx(0.6717, abs=1e-4)
    assert r["m_star"]["loop"] == pytest.approx(0.5740, abs=1e-4)
    assert r["band_pass"] == pytest.approx({"heat kernel": 0.0661, "exponential in |D|": 0.1009,
                                            "Gaussian in D^2": 0.0558, "eigenvalue counting": 0.0522,
                                            "fermionic entropy": 0.0726}, abs=1e-3)
    assert r["anchor_edge"] == pytest.approx(0.999945, abs=2e-6) and r["labelled_at_edge"] == pytest.approx(-0.174, abs=2e-3)
    assert r["kernel_sigma"] == pytest.approx(25.8333, abs=1e-3) and r["kernel_shape"] == pytest.approx(-2.0, abs=1e-3)
    assert r["static_sign"] == pytest.approx(1.0, abs=1e-3)
    assert r["zero"] is None and r["kernel_limit"] == 65 / 8
    assert r["mass_response_one"] == 0.0 and r["mass_response_massive"] == pytest.approx(-0.988, abs=1e-3)
    assert r["regulator"] == pytest.approx(0.4932, abs=1e-3)
    assert r["q_departure"] == pytest.approx([0.130, 0.432], abs=2e-3)
    assert r["n_edge"] == pytest.approx(139.31, abs=0.05) and r["zero_free"]
    assert r["q_zero_label"] == pytest.approx(1.295, abs=2e-3)
    assert r["xi_band"]["no loop"] == pytest.approx((-0.1781, 0.1074), abs=1e-4)
    assert r["xi_band"]["loop"] == pytest.approx((-0.1928, 0.0928), abs=1e-4)
    assert r["xi_edge"] == pytest.approx(0.01508, abs=1e-4)
    assert [rd for rd, v in r["readings"].items() if all(v)] == ["amount"]
    assert r["share_tolerance"] == pytest.approx(5.3e-6, abs=2e-7) and r["one_seven"] == pytest.approx(2.0, abs=0.1)
    assert r["share_quantum"] / r["share_tolerance"] > 3e4
    assert r["readings"]["light"] == (False, True, True) and r["readings"]["boundary"] == (True, False, False)
