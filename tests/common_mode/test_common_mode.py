"""Freezes the common mode's node (interference/root.py, THE COMMON MODE): G M_b^2 = 3 Q_l - 1.

The node has no base and two depth four edges, the powers of the
relative channels q = 1, 2 of the lepton circulant on its depth four
rates, over the q = 0 power.  Each is 1/2 on the searchlight's
circulant.  On the depth four web they sum to
P v^2 P = 3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2
(tests/probes/common_mode_count.py).  G reads the ruler through it.
The node carries no weight of its own, and the count is fixed with no
measured G.
"""

import math

import pytest

G_REC, G_UNC = 6.67430e-11, 0.00015e-11              # CODATA recommended value and its uncertainty


def _pull(G):
    return (G - G_REC) / G_UNC


def test_the_node_is_base_zero_and_two_depth_four_counts():
    """The node's grammar: no base, two depth four vents with parent
    eq_loop, each 1/2 on the searchlight's circulant, web depth four.
    The ancestry rule passes, and its value is 3 Q_l - 1 = 1 - 7.50 L."""
    from interference import root as r

    led = r.WEB["common_mode"]
    assert led.base == 0.0 and led.mode == "add" and led.web_depth == 4
    assert len(led.terms) == 2 and all(t.depth == 4 and t.kind == "vent" and t.parent == "eq_loop"
                                       for t in led.terms)
    searchlight = dict(r.WEB.state, eq_loop=0.0)
    assert [t.value(searchlight) for t in led.terms] == pytest.approx([0.5, 0.5], abs=1e-14)
    r.WEB.validate_depth_ancestry()
    L = r.WEB.state["eq_loop"]
    rates = [g * g * (1 + r.lepton_carrier_term(r.WEB.state, g)) for g in r.LEPTON_GAPS.values()]
    Q = sum(rates) / sum(math.sqrt(x) for x in rates) ** 2
    assert r.COMMON_MODE == r.WEB.state["common_mode"] == pytest.approx(3 * Q - 1, abs=1e-14)
    assert (r.COMMON_MODE - 1) / L == pytest.approx(-7.5022, abs=1e-4)
    assert r.COMMON_MODE - 1 == pytest.approx(-1.011956e-5, rel=1e-5)


def test_G_and_its_pull():
    """Through CODATA's printed Planck mass G = 6.673936e-11, -2.43
    sigma of the recommended value.  The ruler does not move: M_Pl is
    the anchor inversion, +22.23 ppm above CODATA's Planck mass, and the
    Planck mass G reads, sqrt(hbar c/G), sits 27.29 ppm above it."""
    from interference import root as r

    assert r.G_PRED == pytest.approx(6.673935680964234e-11, abs=1e-22)
    assert r.G_PRED == pytest.approx(r.COMMON_MODE * r.G_DEPTH_ONE, rel=1e-15)
    assert _pull(r.G_PRED) == pytest.approx(-2.429, abs=1e-3)
    assert 1e6 * (r.M_Pl_GeV / r.M_PL_CODATA_GEV - 1) == pytest.approx(22.23, abs=0.01)
    assert 1e6 * (r.M_Pl_GeV / math.sqrt(r.COMMON_MODE) / r.M_PL_CODATA_GEV - 1) == pytest.approx(27.29, abs=0.01)
    assert r.RHO_CRIT_MAPPED == pytest.approx(3 * r.H_0_GEV**2 * r.M_Pl_GeV**2 / (8 * math.pi) / r.COMMON_MODE,
                                              rel=1e-15)


def test_the_node_moves_no_other_node(res):
    """Switch the edges off and on, on the solved web of the full
    chain.  Every node except the common mode keeps its value to the
    last bit.  With its edges off the common mode falls to its base, 0.
    No web node reads it back."""
    from interference import root as r

    web, led = r.WEB, r.WEB["common_mode"]
    on = dict(web.state)
    edges = list(led.terms)
    try:
        led.terms.clear()
        off = dict(web.solve())
    finally:
        led.terms.extend(edges)
        web.solve()
    assert off["common_mode"] == 0.0 and web.state["common_mode"] == on["common_mode"]
    moved = [k for k in on if k != "common_mode" and off[k] != on[k]]
    assert moved == [] and len(on) >= 28


def test_the_alternatives_and_what_the_measurement_resolves():
    """Against the recommended value the count gives -2.43.  A vertex
    weight w applied to the depth one whole gives -1.92 for the rule's
    1, -0.36 for the knot's traversal 27, -1.86 for the orientation
    count 2, -3.60 for the amplitude's sign -27 and -1.98 for no edge.
    Read as a conservation rule ("the common mode gains what closures
    lose"), the depth two factor (1 - alpha/2pi) would put G at +49.7.
    The recommended G needs a second order slope of 33.0 +- 16.7, and
    the count's -7.50 lies 2.43 of that width below it."""
    from interference import root as r

    L = r.WEB.state["eq_loop"]
    pulls = {w: _pull(r.G_DEPTH_ONE * (1 + w * L)) for w in (1, 27, 2, -27, 0)}
    assert [round(pulls[w], 2) for w in (1, 27, 2, -27, 0)] == [-1.92, -0.36, -1.86, -3.60, -1.98]
    assert _pull(r.G_DEPTH_ONE * (1 + math.sqrt(L))) == pytest.approx(49.7, abs=0.1)
    need = (G_REC / r.G_DEPTH_ONE - 1) / L
    width = G_UNC / r.G_DEPTH_ONE / L
    slope = (r.COMMON_MODE - 1) / L
    assert need == pytest.approx(32.97, abs=0.01) and width == pytest.approx(16.66, abs=0.01)
    assert (slope - need) / width == pytest.approx(-2.43, abs=0.01)
    assert _pull(r.G_PRED) == pytest.approx(_pull(r.G_DEPTH_ONE * (1 + slope * L)), abs=1e-3)
