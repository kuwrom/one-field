"""The vertex composition rule read on the common mode, with one classification, and why it is not the node's reading.

THE QUESTION. Every depth 4 reading of the closed carrier L = (alpha/2pi)^2
has, by the vertex composition rule, "at each endpoint, THE SAME VERTEX
WEIGHT that node's channel already carries at depth 1, echoes compose
existing vertices, they never invent new ones" (root.py). The rule
gives the coupling 2, the quark action 16 and the lepton amplitudes
27/2. What does the same rule give the common mode, over the whole
frozen vertex menu, with no measured G? Is that weight the node's
reading?

PREMISES. P1 the vertex composition rule and its frozen menu, e-side
{1, 2, 8/3, 1/2} by q-side {1, 6, 27, 26, 9, 3/5, 182}
(tests/probes/wiring_theorem_probe.py PART D). P2 the common mode's
edges: the powers of the two relative channels, 1/2 each on the
searchlight's circulant, 1 when read whole. It carries no charge vertex
and no F4-side vertex (root.py, THE COMMON MODE). P3 the common mode is
the invariant component. On the Z3 side that is P. At the F4 end it is
the F4-invariant part of the Albert traversal, J3(O) = 1 + 26 under F4
(protected_forgetting.py; wiring_theorem_probe.py, the unital transport
T = P1 + P26). P4 a knot carries names and reads the whole traversal.
m_mu/m_e holds the lepton weight at 13.50 +- 0.03
(frame_coefficients.py F3).

STATEMENTS (the numbers below check each one).

  W1 (the e-side). The common mode's depth one channel carries its two
     vents, read whole: 1. The charge trace 8/3 needs a charge vertex,
     and the neutral node has none. The orientation count 2 is the
     rule's weight for a node with no vertices, and the common mode has
     two. W(fund) = 1/2 is one channel's share, and the common mode
     reads the whole (boundary_lattice.py). One entry of four stands: 1.
  W2 (the q-side). The common mode's reading of the unital Albert
     traversal is its invariant part, P1 T P1, trace 1. The whole trace
     27 includes the 26, which F4 moves, and is a knot's reading (P4).
     The Casimir 6, the dual Coxeter 9 and the weight 3/5 are vertices
     of other ledgers. 26 is the moved part alone. 182 is the amount's
     count. One entry of seven stands: 1.
  W3 (the weight and its sign). Of the 28 products, one reading stands:
     1 x 1 = 1. The node's edges are powers and the weight a count, so
     the rule's term is +L: the rule reads the common mode as 1 + L.
  W4 (the test). With weight 1, G = 6.674012e-11, -1.92 standard
     deviations from CODATA's recommended value. The weight 27, the knot's
     traversal read by the common mode, gives -0.36. The coupling's
     orientation count 2 gives -1.86. The amplitude's sign, -27, gives
     -3.60. The measurement asks for 33 +- 17 and does not resolve
     weights between 0 and 2.
  W5 (the node's reading). The rule composes a vertex weight onto the
     response a node's channel already carries. The common mode has no
     response of its own. Its edges are the relative powers of the
     lepton circulant on the depth four rates. The carrier reaches them
     through the three lepton terms, -7.50 L. The node's reading is that
     count, G M_b^2 = 3 Q_l - 1. A vertex weight added beside it would
     read the carrier a second time (tests/probes/common_mode_count.py
     C2, C3: K6 and K7 fail).

RESULT. Read with one classification, the vertex composition rule gives
the common mode the weight 1: its own vents read whole, and the
invariant part of the Albert traversal. The 27 would need the common
mode to read the 26, a name, as a knot does. The weight 1 is the
rule's, not the node's reading. The node has no response of its own
for a vertex weight to compose onto. It reads the carrier once,
through the lepton circulant it counts
(tests/probes/common_mode_count.py). The rule's weight 1 would put G
at -1.92 standard deviations from the recommended value. The count
puts it at -2.43. Both distances depend on the rule used to combine
the measurements (tests/audit/boundary_reading_scan.py).

Usage: python3 tests/probes/common_mode_weight.py
"""

import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402

G_CODATA = (6.67430e-11, 0.00015e-11)
E_SIDE = {"1": Fraction(1), "2 (orientation)": Fraction(2), "8/3 (charge trace)": Fraction(8, 3),
          "1/2 (W(fund))": Fraction(1, 2)}
Q_SIDE = {"1": Fraction(1), "6 (C2(26))": Fraction(6), "27 (dim J3(O))": Fraction(27),
          "26 (n26)": Fraction(26), "9 (h_dual F4)": Fraction(9), "3/5 (h26)": Fraction(3, 5),
          "182 (N_bridge)": Fraction(182)}


def pull(weight):
    G = r.G_DEPTH_ONE * (1 + weight * r.WEB.state["eq_loop"])
    return G, (G - G_CODATA[0]) / G_CODATA[1]


def run(report=print):
    out = {}
    assert len(E_SIDE) * len(Q_SIDE) == 28

    # W1: the e-side entries the common mode's own channel carries, on the searchlight's circulant.
    searchlight = dict(r.WEB.state, eq_loop=0.0)
    vents = [round(t.value(searchlight), 15) for t in r.WEB["common_mode"].terms]
    whole = sum(vents)
    e_ok = {k: v for k, v in E_SIDE.items() if v == Fraction(whole).limit_denominator()}
    assert vents == [0.5, 0.5] and list(e_ok) == ["1"]

    # W2: the invariant part of the unital Albert traversal.
    n = int(r.n26) + 1
    p1 = np.zeros((n, n))
    p1[0, 0] = 1.0
    T = p1 + (np.eye(n) - p1)
    invariant = float(np.trace(p1 @ T @ p1))
    q_ok = {k: v for k, v in Q_SIDE.items() if v == Fraction(round(invariant))}
    assert invariant == 1.0 and np.trace(T) == 27 and list(q_ok) == ["1"]

    # W3: the reading that stands.
    stands = [(e, q) for e in e_ok for q in q_ok]
    weight = float(E_SIDE[stands[0][0]] * Q_SIDE[stands[0][1]])
    assert stands == [("1", "1")] and weight == 1.0
    out.update(vents=vents, invariant=invariant, stands=stands, weight=weight)
    report(f"W1-W3 e-side {list(e_ok)}, q-side {list(q_ok)} (P1 T P1 = {invariant:.0f} of trace {np.trace(T):.0f}): "
           f"weight {weight:.0f}, the rule's reading 1 + L")

    # W4: the test.
    readings = {"rule, 1": 1.0, "knot traversal, 27": 27.0, "orientation, 2": 2.0, "amplitude sign, -27": -27.0,
                "no edge, 0": 0.0}
    tests = {k: pull(w) for k, w in readings.items()}
    assert abs(tests["rule, 1"][1] + 1.92) < 0.01 and abs(tests["knot traversal, 27"][1] + 0.36) < 0.01
    assert abs(tests["orientation, 2"][1] + 1.86) < 0.01 and abs(tests["amplitude sign, -27"][1] + 3.60) < 0.01
    out.update(tests=tests)
    report("W4 " + "; ".join(f"{k}: G = {g:.6e} ({p:+.2f} sigma)" for k, (g, p) in tests.items()))

    # W5: the node's reading is the count, which the rule's weight does not give.
    L = r.WEB.state["eq_loop"]
    node = r.WEB["common_mode"]
    slope = (r.COMMON_MODE - 1.0) / L
    assert node.base == 0.0 and all(t.parent == "eq_loop" and t.depth == 4 for t in node.terms)
    assert abs(slope + 7.502) < 0.001 and abs(r.COMMON_MODE - (1 + weight * L)) > 8 * L
    assert abs(r.G_PRED - tests["rule, 1"][0]) > 5e-16
    out.update(node_slope=slope)
    report(f"W5 the node's count moves by {slope:+.3f} L through the lepton terms. The rule's weight reads +{weight:.0f} L")
    return out


if __name__ == "__main__":
    run()
