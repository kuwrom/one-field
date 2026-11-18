"""Freezes tests/probes/weak_winding.py: what the weak winding does to the baryon and lepton labels."""

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import weak_winding
    return weak_winding.run(report=lambda *a: None)


def test_the_winding_moves_both_labels_by_three_and_is_neutral(r):
    """W1, W2: one unit changes B and L by d11 = 3 each and the charge
    by 0. On the menu d10 = 1 to 4 it is neutral for d10 = 2 alone."""
    assert r["step"] == {"B": 3, "L": 3, "Q": 0}
    assert r["neutral_d10"] == [2]


def test_the_proton_and_the_neutrino_gate_are_untouched(r):
    """W3 to W5: a single proton cannot decay through the winding, the
    rate factor is 10^-161, and the winding gives no Majorana term."""
    assert r["proton_nearest"] == -2 and r["lightest_system"] == 2
    assert r["alpha_W"] == pytest.approx(0.0339, abs=2e-4)
    assert r["log10_rate"] == pytest.approx(-161.0, abs=0.5)
    assert r["majorana_units"] == [] and r["label_gate_same"]


def test_the_sum_of_b_minus_l_over_the_content(r):
    """W6: no edge and no weak winding changes B - L. Summed over the
    left handed content it is -1 with two right handed neutrinos."""
    assert r["b_minus_l_sum"] == -1
