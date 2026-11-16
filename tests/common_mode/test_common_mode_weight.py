"""Freezes tests/probes/common_mode_weight.py: the vertex composition rule applied to the common mode, and why its weight is not the node's value."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def out():
    from probes import common_mode_weight
    return common_mode_weight.run(report=_silent)


def test_w5_the_nodes_reading_is_the_count(out):
    """The node has no response of its own for a vertex weight to act
    on.  It has base 0 and two depth four edges on the lepton circulant.
    The carrier moves its count by -7.50 L through the lepton terms, not
    by the rule's +1 L.  G_PRED comes from the count, not from the
    rule."""
    from interference import root as r

    node = r.WEB["common_mode"]
    assert node.base == 0.0 and [t.depth for t in node.terms] == [4, 4]
    assert out["node_slope"] == pytest.approx(-7.502, abs=1e-3)
    assert r.G_PRED == pytest.approx(6.673936e-11, abs=1e-17)
    assert abs(r.G_PRED - out["tests"]["rule, 1"][0]) > 5e-16
