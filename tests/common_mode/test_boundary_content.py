"""Freezes tests/probes/boundary_content.py: the labelled bridge on one list read by the running and the stiffness."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_boundary_content_the_labelled_bridge_puts_gate_s_on_the_regulator_shape():
    """Read with its searchlight label on one list, the bridge and the 7
    G2 scalars cannot be light (b0 = 19/3, 21/2: the anchor moves by
    1.6e13 and 2.02).  Without them the covariant contact stiffness is
    -4.17.  The counted ledger, with the bridge read as amount, gives 27.33
    and puts the anchored regulator at 0.48 M_b.  The labelled bridge at
    the boundary needs 0.359 of a massless row's weight for an
    attractive stiffness whose anchored regulator lies inside the
    boundary (0.143 for the sign alone).  The sharp cutoff gives 0.148
    (Sigma <= +0.158, short of 2 pi), and exp(-t0/t) gives 0.280
    (attractive). The G2-singlet content the
    repository names adds 2.01 of the 6.12 at its most favourable
    weights."""
    from probes import boundary_content as bc

    r = bc.run(report=_silent)
    assert r["index_E8"] == 30
    assert r["sigma_booked"] == pytest.approx(82 / 3, abs=1e-12)
    assert r["v_break"] == pytest.approx(1.43e-17, rel=0.01)
    assert r["light"]["bridge"]["M_Pl"] == pytest.approx(1.58e13, rel=0.01)
    assert r["light"]["G2 scalars"]["M_Pl"] == pytest.approx(2.024, abs=0.002)
    assert r["light"]["G2 scalars"]["M_Pl_mu_star"] == pytest.approx(1.823, abs=0.002)
    assert r["sigma_absent"][1.0] == pytest.approx(-25 / 6, abs=1e-12)
    assert r["sigma_absent"][0.0] == pytest.approx(83 / 6, abs=1e-12)
    assert r["regulator_booked"] == pytest.approx(0.4795, abs=1e-3)
    assert r["anchor_edge"] == pytest.approx(0.999945, abs=2e-6)
    assert r["bridge_cap"] == pytest.approx(4.325, abs=1e-3)
    assert r["sigma_inside"] == pytest.approx(0.158, abs=1e-3)
    assert r["needed"] == pytest.approx(6.125, abs=2e-3) and r["needed_scalars"] == 37
    assert r["needed_loop"] == pytest.approx(8.791, abs=2e-3)
    assert r["w_need"] == pytest.approx(0.3588, abs=1e-3) and r["w_need_loop"] == pytest.approx(0.4503, abs=1e-3)
    assert r["w_sharp"] == pytest.approx(0.1485, abs=1e-4) and r["w_smooth"] == pytest.approx(0.2797, abs=1e-4)
    assert r["sigma_smooth"] == pytest.approx(3.981, abs=2e-3)
    assert r["w_sign"] == pytest.approx(0.1431, abs=1e-3)
    assert r["stack"] == pytest.approx(2.007, abs=2e-3)
