"""Freezes tests/probes/bridge_b0_check.py: the bridge's O(g^2) piece is the
G2 field's own Maxwell term, so the amount runs nothing and b0 = 32/3 holds."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_bridge_b0_the_O_g2_piece_is_the_gauge_row_not_a_portal():
    """The O(F^2) density of the Hermitian DBI form at the empty field
    vacuum h = c eta is the G2 field's own Maxwell term: the c^2 from the
    determinant and the 1/c^2 from index raising cancel.  So the piece
    does not depend on c or U to order F^2 (determinant_vectors V5).  It
    is already counted as the G2 gauge row, not a new matter gauge
    portal.  The amount's running share is therefore 0, and b0 = 32/3 on
    the one list.  If the bridge is read as light, the label share -13/3
    gives b0 = 19/3 and M_Pl moves by x1.58e13, so the anchor fails.  A
    surviving portal would take that branch, and this test closes it."""
    from probes import bridge_b0_check

    r = bridge_b0_check.run(report=_silent)

    # B1 the analytic cancel: L_F2 = -1/4 F^2 for every c.
    assert all(v == pytest.approx(-0.25, abs=1e-14) for v in r["cancel"].values())

    # B2 c- and U-independence of the second F derivative.
    assert r["worst_c"] < 1e-5 and r["worst_u"] < 1e-5

    # B3 b0 on the one list: amount share 0.
    assert r["b0_amount"] == "32/3" and r["share"] == "0"

    # B4 the counterfactual: the light reading and its kill.
    assert r["b0_light"] == "19/3"
    assert 1e13 < r["M_Pl_move"] < 2e13
