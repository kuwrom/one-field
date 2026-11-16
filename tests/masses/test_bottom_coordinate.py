"""Freezes tests/probes/bottom_coordinate.py: the bottom mass's relation and coordinate."""

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import bottom_coordinate
    return bottom_coordinate.run(report=lambda *a: None)


def test_the_lever_and_the_relation_test(r):
    """39 MeV per GeV of top. The relation on measured inputs is
    4197.8 +- 14.6 MeV, +0.91 sigma of PDG."""
    assert r["per_gev"] == pytest.approx(39.1, abs=0.3)
    assert r["fed"] == pytest.approx(4197.8, abs=0.1) and r["fed_pull"] == pytest.approx(0.91, abs=0.01)


def test_the_references_bracket_the_engine(r):
    p = r["pulls"]
    assert p["PDG m_b(m_b)"] == pytest.approx(1.54, abs=0.01)
    assert p["FLAG Nf=2+1+1"] == pytest.approx(-0.44, abs=0.01)
    assert p["FLAG Nf=2+1"] == pytest.approx(1.14, abs=0.01)
