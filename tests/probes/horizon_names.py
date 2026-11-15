"""What a name is, and why the bridge carries none: the premise T1 in two parts.

THE QUESTION. The lane's premise T1 reads "The common mode past the
horizon carries no names" (notes/open_derivation_tasks.md), and six
readings rest on it: the boundary reading c, the bridge's zero share of
b0, the amount's scalar weight, C_EFF, minimal coupling and the dust
reading. The repository uses "name" in several senses. What single
sense do the readings need, and does the bridge's half follow from the
algebra?

PREMISES. P1 a record is a finite dimensional unital real composition
algebra, and records are primitive (interference/__init__.py, THE
ONTOLOGY). P2 the conformal embedding E8(1) > G2(1) x F4(1) with
248 = (14,1) + (1,52) + (7,26) (root.py). P3 the E8 root lattice
(interference/wiring.py).

THE DEFINITION. A name is a label of the searchlight's data that tells
closures apart: which representation or channel of the algebra a state
sits in. An amount is a number read off the whole: a count, a weight, a
total. A reading reads no name when it reads no representation or
channel label, and it may read amounts.

STATEMENTS (H3 is computed below; H1 and H2 read the cited probes).

  H1 (the six readings in that sense). The boundary reading reads the
     whole set of closures, kept by every symmetry
     (boundary_lattice.py L4, L5). A share of b0 is a Dynkin index.  It
     tells G2 representations apart, so it is a name, and the amount has
     no share (gravity_premises.py Q1). The scalar weight is the reading kept by
     relabelling the components (unnamed_metric.py). C_EFF reads the
     whole second moment (frame_coefficients.py F1). The amount's
     coupling involves no vertex, and every vertex carries a name
     (boundary.py, THE AMOUNT'S COUPLING). The dust reading needs a
     momentum mode to count as a name (T10).  It is the only reading the
     definition does not settle, so it keeps its own premise.  A
     momentum mode labels a representation of translations, not of the
     searchlight's data.
  H2 (gate U in that sense). The weight and the count of a lump are
     both amounts. Gate U compares two readings of amounts: a bias
     that reads the weight keeps the ratio of two clocks, and one that
     reads the count moves it by 0.297 of itself
     (gates_as_records.py U). The count is not a name.  U's failing
     case uses the wrong amount.
  H3 (the bridge is a multiplicity of the unit). The E8 lattice is even
     and unimodular, Gram determinant 1, so E8(1) has one primary, the
     vacuum. Its weight one states are the 248 currents, and the bridge's
     7 x 26 = 182 sit among them beside 14 and 52. Read in the ambient
     theory no current has a label of its own, the bridge among them.
     Read by G2 x F4, the searchlight's split, it is the bifundamental,
     a name. A label on the bridge as one object needs a record holding a
     joint (7, 26) label, and O tensored with a nontrivial algebra has
     zero divisors, so no composition algebra holds it (P1).

RESULT. One definition serves five of the six readings: a name is a
representation or channel label of the searchlight's data, and amounts
are numbers read off the whole. The bridge
carries no name because the ambient theory has one primary and no
composition algebra holds a joint label. The dust reading needs momentum
modes to be names, which stays the premise T10.

Usage: python3 tests/probes/horizon_names.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402
from interference import wiring  # noqa: E402


def run(report=print):
    out = {}
    G = wiring._gram_E8()
    det = round(float(np.linalg.det(G)))
    even = all(int(G[i, i]) % 2 == 0 for i in range(8))
    roots = wiring.e8_roots()
    weight_one = len(roots) + 8
    split = (r.dim_G2, r.dim_F4, r.N_bridge)
    assert det == 1 and even and weight_one == 248 == sum(split) == r.dim_E8
    J = np.array([[0.0, -1.0], [1.0, 0.0]])            # an imaginary unit of O, commuting with i
    zero = (np.eye(2) + 1j * J) @ (np.eye(2) - 1j * J)
    assert np.allclose(zero, 0) and not np.allclose(np.eye(2) + 1j * J, 0)
    out.update(det=det, even=even, weight_one=weight_one, split=split, zero_divisor=True)
    report(f"H3 E8 Gram determinant {det}, even {even}: one primary. Weight one states {weight_one} = "
           f"{split[0]} + {split[1]} + {split[2]}; (1 + i e)(1 - i e) = 0 in O tensored with C")
    return out


if __name__ == "__main__":
    run()
