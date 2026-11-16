"""The G-path lane exhaustion (tests/audit/g_lane_exhaustion.py).

What is frozen: across the full product of the G-path menus
(4,305,920 lanes, with c2 over the rival readings of
frame_coefficients.py F2):
  * G alone closes many lanes;
  * G_F and m_mu fix C_EFF, c2, w4 and w16 without G;
  * exactly one lane closes on G, G_F and m_mu together: the canonical
    one.  The next lane that closes on G_F and m_mu misses G by
    988 sigma.
"""

import os
import sys

import pytest

_here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_here, "audit"))
sys.path.insert(0, os.path.join(_here, "probes"))
import g_lane_exhaustion as gle  # noqa: E402


@pytest.fixture(scope="module")
def ex():
    return gle.run(report=lambda *a, **k: None)


def test_g_alone_does_not_select(ex):
    assert ex["n_lanes"] == 4 * 4 * 2 * 4 * 8 * 5 * 29 * 29 == 4305920
    assert ex["close_G"] > 100


def test_other_observables_pin_the_fine_coefficients(ex):
    assert ex["pinned"] == ({"1/2"}, {"h_bridge"}, {"27/2"}, {("16",)})


def test_exactly_one_lane_closes(ex):
    assert len(ex["close_all"]) == 1
    lane = ex["close_all"][0]
    assert lane[:7] == ("1/2", "32/3", "P h/2pi", "action", "1", "h_bridge", "27/2")
    assert lane[-1] == pytest.approx(-2.43, abs=0.01)
    assert ex["runner_sigma"] == pytest.approx(988, abs=1)   # exact pull
