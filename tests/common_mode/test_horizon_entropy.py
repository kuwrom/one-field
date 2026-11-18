"""Freezes tests/probes/horizon_entropy.py: a horizon's entropy on the one list against its rows' entanglement."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_horizon_entropy_is_the_rows_entanglement_plus_the_vector_edge_term():
    """The cone's image sums continue to an entropy weight of 1/6 per
    scalar and 1/12 per Weyl state, which are their heat kernel weights.
    The one list's entanglement is 44.1667, equal to its stiffness at
    contact 0.  Across the contacts this is 1.688, 1 and 0.710 of the
    horizon's A/(4G).  At the covariant contact the difference is the
    vectors' -18.  The amount's 30.3333 is all entanglement.  The
    label's xi would add 1.207 that no image sum contains.  The
    nameable rows alone have entanglement 13.83 against a stiffness of
    -4.17."""
    from probes import horizon_entropy as he

    r = he.run(report=_silent)
    assert r["scalar_weight"] == pytest.approx(1 / 6, abs=1e-10)
    assert r["spinor_weight"] == pytest.approx(1 / 12, abs=1e-10)
    assert r["entanglement"] == pytest.approx(44.1667, abs=1e-4) and r["contact_term"] == pytest.approx(-18.0)
    assert [r["share"][c] for c in (1.0, 0.0, -1.0)] == pytest.approx([1.6879, 1.0, 0.7105], abs=1e-3)
    assert r["amount"] == pytest.approx(30.3333, abs=1e-4) and r["label_gap"] == pytest.approx(1.2069, abs=1e-3)
    assert r["named"] == pytest.approx((83 / 6, -25 / 6), abs=1e-12)
