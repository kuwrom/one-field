"""Freezes tests/probes/horizon_amount.py: gravity's sign and size read only the amount past the searchlight's horizon."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight():
    """The nameable rows give -4.17 at the covariant contact (-6.83 with
    the metric loop), which is repulsive.  Neither an octonion's worth of
    real scalars nor the named G2 singlets make them attractive.  A light
    amount of 78.7 to 390.5 real scalar units passes S and M at every
    contact with the loop (62.7 to 374.5 without), for any cutoff shape
    and placement.  The searchlight's 182 lies inside this window with a
    factor of about two to spare.  The label's coupling on the amount
    moves the edges by 4 percent.  With a mass, the amount keeps both
    gates below 0.685 M_b (0.586 M_b with the loop) under the least
    favourable monotone regulator.  Band pass shapes at the spectral
    amplitude need 0.05 to 0.10 M_b.  The labelled bridge keeps the
    anchor only above 0.999945 M_b, where the ledger's cutoff leaves it
    at +0.159, below the 2 pi that S needs.  The one list's kernel is
    attractive and keeps the Einstein sign on the axis, with the limit
    63/8.  A zero falls inside
    the boundary only below 140.1 units.  The label's coupling on the
    amount would put one at 1.268 M_b.  At the count the gates hold for
    -0.176 < xi < 0.109, and the zero stays outside the boundary for
    xi < 0.0147.  Of the five readings the list can give the bridge,
    only the amount keeps b0 and passes S and M."""
    from probes import boundary_content as bc
    from probes import horizon_amount as ha

    r = ha.run(report=_silent)
    assert ha.NAMED_SINGLETS == pytest.approx(bc.run(report=_silent)["stack"], abs=2e-3)
    assert r["named"] == pytest.approx({1.0: -25 / 6, 0.0: 83 / 6, -1.0: 191 / 6}, abs=1e-12)
    assert r["named_loop"] == pytest.approx(-6.8333, abs=1e-4) and r["octonion"] == pytest.approx(4 / 3)
    assert r["attract"] == pytest.approx(25.0, abs=1e-2) and r["attract_loop"] == pytest.approx(41.0, abs=1e-2)
    win = r["window"]
    assert win["no loop"] == pytest.approx((62.70, 374.49), abs=1e-2)
    assert win["loop"] == pytest.approx((78.70, 390.49), abs=1e-2)
    assert win["no loop, label's coupling"] == pytest.approx((65.30, 390.00), abs=1e-2)
    assert win["loop, label's coupling"] == pytest.approx((81.96, 406.67), abs=1e-2)
    assert r["margins"]["loop"] == pytest.approx((2.313, 2.146), abs=1e-3)
    assert r["margins"]["no loop"] == pytest.approx((2.903, 2.058), abs=1e-3)
    assert list(r["sigma"].values()) == pytest.approx([26.1667, 44.1667, 62.1667], abs=1e-4)
    assert r["t0"][1.0] == pytest.approx(4.1646, abs=1e-4)
    assert r["m_star"]["no loop"] == pytest.approx(0.6850, abs=1e-4)
    assert r["m_star"]["loop"] == pytest.approx(0.5856, abs=1e-4)
    assert r["band_pass"] == pytest.approx({"heat kernel": 0.0671, "exponential in |D|": 0.1030,
                                            "Gaussian in D^2": 0.0564, "eigenvalue counting": 0.0526,
                                            "fermionic entropy": 0.0738}, abs=1e-3)
    assert r["anchor_edge"] == pytest.approx(0.999945, abs=2e-6) and r["labelled_at_edge"] == pytest.approx(0.159, abs=2e-3)
    assert r["kernel_sigma"] == pytest.approx(26.1667, abs=1e-3) and r["kernel_shape"] == pytest.approx(-2.0, abs=1e-3)
    assert r["static_sign"] == pytest.approx(1.0, abs=1e-3)
    assert r["zero"] is None and r["kernel_limit"] == 63 / 8
    assert r["mass_response_one"] == 0.0 and r["mass_response_massive"] == pytest.approx(-0.979, abs=1e-3)
    assert r["regulator"] == pytest.approx(0.4900, abs=1e-3)
    assert r["q_departure"] == pytest.approx([0.129, 0.430], abs=2e-3)
    assert r["n_edge"] == pytest.approx(140.09, abs=0.05) and r["zero_free"]
    assert r["q_zero_label"] == pytest.approx(1.268, abs=2e-3)
    assert r["xi_band"]["no loop"] == pytest.approx((-0.1763, 0.1092), abs=1e-4)
    assert r["xi_band"]["loop"] == pytest.approx((-0.1909, 0.0946), abs=1e-4)
    assert r["xi_edge"] == pytest.approx(0.01469, abs=1e-4)
    assert [rd for rd, v in r["readings"].items() if all(v)] == ["amount"]
    assert r["share_tolerance"] == pytest.approx(5.3e-6, abs=2e-7) and r["one_seven"] == pytest.approx(2.0, abs=0.1)
    assert r["share_quantum"] / r["share_tolerance"] > 3e4
    assert r["readings"]["light"] == (False, True, True) and r["readings"]["boundary"] == (True, False, True)
