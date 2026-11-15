"""Mechanism probes of the foundations branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##

def test_wiring_found_not_stated():
    """E8 > G2(1) x F4(1) at k = 3 is the unique zero residual wiring.
    The dictionary is read off the winning wiring. It is never an input."""
    import wiring
    r = wiring.run(report=_silent)
    assert r["lattice"] == "E8"
    assert r["pair"] == "G2(1) x F4(1)"
    assert r["k"] == 3
    assert (r["n7"], r["n26"]) == (7, 26)
    assert abs(r["d10"] - 2.0) < 1e-9 and abs(r["d11"] - 3.0) < 1e-9


def test_normalized_current_response_wired():
    """The selected h=1 simple current and its S-matrix channel weight
    generate alpha_G2(M_Pl). root.py uses the same wiring record."""
    import wiring
    import root

    w = wiring.canonical()
    assert w["k"] == 3 and w["h_J"] == 1.0
    assert w["pointed_orbit"] == ((0, 0), (3, 0), (0, 3))
    assert abs(w["D2_tot"] - 36.0) < 1e-12
    assert abs(w["pointed_weight"] - 1.0 / 12.0) < 1e-12
    assert abs(w["current_response"] - 1.0 / (2.0 * math.pi)) < 1e-12
    assert w["current_response"] == wiring.normalized_current_response(
        root.h_bridge)
    with pytest.raises(ValueError):
        wiring.normalized_current_response(2.0)
    assert abs(w["alpha_G2_Pl"] - 1.0 / (24.0 * math.pi)) < 1e-12
    assert w["N_bridge"] == 7 * 26 == root.N_bridge
    assert abs(w["bridge_mode_moment"] - 0.5) < 1e-12
    assert w["bridge_mode_moment"] == wiring.bridge_mode_second_moment(
        math.sqrt(w["d10"]), w["d10"] / w["d11"]**2,
        w["n7"], w["n26"])
    assert float(root.E_v2) == pytest.approx(
        w["bridge_mode_moment"], abs=1e-12)
    assert root.WIRING_DATA == w
    assert root.alpha_G2_Pl == w["alpha_G2_Pl"]


def test_e8_vertex_charge_trace_wired():
    """The one orientation E8 current index is 30 for every unit q."""
    import wiring
    import root
    form = wiring.e8_positive_root_current_form()
    assert form == pytest.approx(30.0 * np.eye(8), abs=1e-12)
    assert len(wiring.e8_roots()) == 240
    assert wiring.e8_vertex_charge_trace() == pytest.approx(30.0)
    diagonal = np.ones(8) / math.sqrt(8.0)
    assert wiring.e8_vertex_charge_trace(diagonal) == pytest.approx(30.0)
    assert root.vertex_charge_trace == root.hv_E8 == root.N_vertex == 30
    assert root.WIRING_DATA["vertex_charge_trace"] == 30


## Statistics: the WZ level by counting, and the winding ##

def test_bookkeeping_level_and_theta():
    """Child anomaly vanishes. WZ level = 3 by counting. The lightest
    Z3 sector changer is h = 2/9 (internal assertions)."""
    from probes import bookkeeping
    bookkeeping.run(report=_silent)


def test_selector_theorem_binary():
    """The idempotent measure space on Z3 is exactly {Haar, delta_e}:
    the confinement dichotomy as a two element theorem."""
    from probes import selector_theorem
    found = selector_theorem.run(report=_silent)
    assert len(found) == 2


def test_conversion_rigidity():
    """Label rigidity (the conversion lemma) closes the even switch
    loophole. Hom(x,y) = delta over all pairs, there is no same
    orientation vertex, and conversion = pair fusion (vent kind). No
    ghosts appear."""
    from probes import conversion_rigidity
    r = conversion_rigidity.run(report=_silent)
    assert r["rigidity"] and r["conversion_dim"] == 1
    assert r["ghost_fence_gev"] == [64, 128, 139, 277]


def test_wiring_theorem_probe():
    """The wiring theorem: Weyl/Kac-Peterson data, the unique (7,26)
    boson, the frozen 4 x 7 response inventory (PART D, the executable
    menu vertex_pinning cites), and the edge classification with zero
    UNCLASSIFIED canon edges."""
    from fractions import Fraction
    from probes import wiring_theorem_probe as wt
    r = wt.run(report=_silent)
    assert r["n_edges"] == 21
    assert r["counts"] == {"EMBEDDING": 12, "INTRA-FACE": 9,
                           "UNCLASSIFIED": 0}
    assert r["albert"]["options"] == {
        0: Fraction(13), 1: Fraction(27, 2)}
    assert r["albert"]["junction"] == 1
    assert r["albert"]["response"] == Fraction(27, 2)


def test_dimensionality_derived():
    """3+1 is DERIVED, not assumed. Three components force a C^2
    relative sector (codim 4, measured: D=3 no zeros, D=4 isolated
    points, D=5 lines). Knots need d_spatial >= 3. Pointlike charge
    events need D = 4. So D = 3+1 is unique, on the gas sign, where the
    relative sector is populated."""
    from probes import dimensionality
    r = dimensionality.run(report=_silent)
    assert r["unique_D"] == 4
    assert r["codim"][4]["cnt_hi"] / r["codim"][4]["cnt_lo"] < 1.5


def test_dimension_by_closure():
    """3+1 read from closure on the committed sign. Circular orbits are
    stable only below four space dimensions. Near circular and eccentric
    orbits close (apsidal angle pi) only at d = 3. A trial state's energy
    is bounded below only below d = 4. The metric free Maxwell term is
    free of U only at D = 4."""
    from probes import dimension_by_closure as dc

    r = dc.run(report=_silent)
    c = r["curvature"]
    assert c[2] > 0 and c[3] > 0 and abs(c[4]) < 1e-6 and c[5] < 0 and c[6] < 0
    for e in (0.01, 0.3, 0.7):
        assert r["apsidal"][(3, e)] == pytest.approx(math.pi, rel=1e-6)
    assert r["apsidal"][(2, 0.01)] == pytest.approx(math.pi / math.sqrt(2), rel=1e-3)
    assert abs(r["apsidal"][(2, 0.7)] - r["apsidal"][(2, 0.01)]) > 0.05
    assert r["bounded"] == {2: True, 3: True, 4: True, 5: False}
    assert r["weyl_exponent"][4] == 0 and all(v != 0 for D, v in r["weyl_exponent"].items() if D != 4)


def test_envelope_universality():
    """The NLS form is DERIVED up to locality. The envelope of a
    microdynamics other than NLS (phi^4 Klein-Gordon) obeys the NLS
    equation to normal form accuracy, with error ~ eps^2."""
    from probes import envelope_universality
    r = envelope_universality.run(report=_silent)
    assert r["ratio"] < 0.45
    assert r["errors"][0.06][0] < 0.025
