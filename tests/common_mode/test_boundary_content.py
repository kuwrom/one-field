"""Freezes tests/probes/boundary_content.py: the labelled bridge on one list read by the running and the stiffness."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_boundary_content_the_labelled_bridge_puts_the_sign_on_the_regulator_shape():
    """Read with its searchlight label on one list, the bridge and the 7
    G2 scalars cannot be light (b0 = 19/3, 21/2: the anchor moves by
    1.6e13 and 2.02).  Without them the covariant contact stiffness is
    -4.5.  The counted ledger, with the bridge read as amount, gives 27
    and puts the anchored regulator at 0.48 M_b.  The labelled bridge at
    the boundary needs 0.370 of a massless row's weight for an
    attractive stiffness whose anchored regulator lies inside the
    boundary (0.154 for the sign alone).  The sharp cutoff gives 0.148
    (Sigma <= -0.175, repulsive), and exp(-t0/t) gives 0.280
    (attractive). The G2-singlet content the
    repository names adds 2.34 of the 6.46 at its most favourable
    weights."""
    from probes import boundary_content as bc

    r = bc.run(report=_silent)
    assert r["index_E8"] == 30
    assert r["sigma_booked"] == pytest.approx(27.0, abs=1e-12)
    assert r["v_break"] == pytest.approx(1.43e-17, rel=0.01)
    assert r["light"]["bridge"]["M_Pl"] == pytest.approx(1.58e13, rel=0.01)
    assert r["light"]["G2 scalars"]["M_Pl"] == pytest.approx(2.024, abs=0.002)
    assert r["light"]["G2 scalars"]["M_Pl_mu_star"] == pytest.approx(1.823, abs=0.002)
    assert r["light"]["G2 scalars"]["eta_B"] == pytest.approx(0.731, abs=0.002)
    assert r["sigma_absent"][1.0] == pytest.approx(-4.5, abs=1e-12)
    assert r["sigma_absent"][0.0] == pytest.approx(13.5, abs=1e-12)
    assert r["regulator_booked"] == pytest.approx(0.4824, abs=1e-3)
    assert r["anchor_edge"] == pytest.approx(0.999945, abs=2e-6)
    assert r["bridge_cap"] == pytest.approx(4.325, abs=1e-3)
    assert r["sigma_inside"] == pytest.approx(-0.175, abs=1e-3)
    assert r["needed"] == pytest.approx(6.458, abs=2e-3) and r["needed_scalars"] == 39
    assert r["needed_loop"] == pytest.approx(9.125, abs=2e-3)
    assert r["w_need"] == pytest.approx(0.3702, abs=1e-3) and r["w_need_loop"] == pytest.approx(0.4618, abs=1e-3)
    assert r["w_sharp"] == pytest.approx(0.1485, abs=1e-4) and r["w_smooth"] == pytest.approx(0.2797, abs=1e-4)
    assert r["sigma_smooth"] == pytest.approx(3.646, abs=2e-3)
    assert r["w_sign"] == pytest.approx(0.1545, abs=1e-3)
    assert r["stack"] == pytest.approx(2.341, abs=2e-3)
