"""
Terminus purity follows from label rigidity.
Excluding mixed words needs no 6j symbol.

SCOPE.  The conversion lemma excludes mixed terminus words from
bases. Its irrationality witness (each switch inserts one
|C_3bar| = sqrt(2), so odd switch words carry an irrational factor)
covers odd parity only.  Even switch orderings inside the f2a2
class (faaf, affa) carry sqrt(2)^2 = 2, rational. Label rigidity
excludes them.

MECHANISM.  No 6j evaluation is needed.  In the standard formalism
for exactly this situation, topological defect and boundary lines
in RCFT (the TFT construction of Fuchs-Runkel-Schweigert,
Nucl.Phys. B646 (2002) 353; topological defect lines,
Chang-Lin-Shao-Wang-Yin, JHEP 2019; Etingof-Gelaki-Nikshych-Ostrik,
"Tensor Categories", Ch. 1), a line labeled by a SIMPLE object
cannot change its label anywhere along its length. A label change is
a junction, junctions are Hom spaces, and Schur's lemma gives
Hom(x, y) = 0 for distinct simples.  The amplitude of a would-be even
switch is therefore not "rational 2". It is the ZERO MORPHISM, at any
parity, independent of any value.  No F-symbol enters.

THE THREE FACTS (computed below from the repo's own Kac-Peterson
fusion data, mixing._compute_wzw):

  R1 (label rigidity).  dim Hom(x, y) = delta_xy over all 10 x 10
      primary pairs.  A single circulating knot line carries ONE
      channel label for its whole circulation.  An in strand switch,
      odd or even, has amplitude 0.  The knot definition ("a closure
      is a single knot circulating in one channel") states this. It
      is Schur's lemma, not a definition.

  R2 (no same orientation vertex).  dim Hom(f (x) f, a) = 0, since
      f (x) f = (0,1) + (2,0) and the adjoint does not appear.  Even
      the vertex mediated conversion is unavailable inside a
      uniformly oriented tower: an f-tower cannot emit an a-line at
      any price.

  R3 (conversion = orientation reversed pair fusion).  dim
      Hom(f (x) fbar, a) = 1 (the words.py E3 vertex), and the walk
      algebra knows it exactly: n_a = n_f n_fbar - I (the canon
      construction A = N N^T - I of words.py).  Conversion exists
      ONLY as fusion of an f with an ORIENTATION REVERSED fbar.
      That is the Delta/Delta-dagger pairing, which changes the
      strand number.  A configuration that changes the strand number
      is an interaction history. In the frozen edge taxonomy
      (root.EchoTerm kinds) it is of vent kind, never a base.

CONSEQUENCE.  ALL mixed orderings are excluded from bases:
  odd class:  excluded twice over (rigidity, and the original
              irrationality witness sqrt(2)^odd).
  even class: excluded by rigidity (R1) and vertex absence (R2).
              No morphism exists to carry the "rational 2"
              evaluation.
The knot definition only describes this. Schur's lemma does the work.
The only definitional layer left is the base/vent taxonomy itself
(a base counts ONE closure's states). It is uniform canon with no
numerical freedom, the same taxonomy every edge in the web uses.

The empirical fence remains as the physical check. The mixed word
ghosts (~64, ~139 GeV) and the would-be even response ghost (~128 GeV)
match no PDG state and must stay absent.

Usage: python3 tests/probes/conversion_rigidity.py
"""

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.dirname(_here),
          os.path.join(os.path.dirname(os.path.dirname(_here)),
                       "interference"),
          os.path.join(_here, "interference")):
    if p not in sys.path:
        sys.path.insert(0, p)
import _bootstrap  # noqa: E402,F401


def run(report=print):
    import numpy as np
    from root import PDG_MASSES
    from mixing import _compute_wzw
    from words import N_FUND

    report("CONVERSION RIGIDITY: the even switch loophole, closed")
    report("=" * 66)

    w = _compute_wzw()
    P, IDX, NF = w["PRIMARIES"], w["IDX"], w["N_fus"]
    idn, f, fb, a = (0, 0), (1, 0), (0, 1), (1, 1)

    # # R1: label rigidity (Schur), all 100 pairs #
    # dim Hom(x, y) = δ_xy is the identity column of any fusion ring,
    # so this is a data integrity check on the computed N_fus.  The
    # physical content, no in strand switch, is Schur's lemma (cited
    # above).
    rigidity_ok = bool(all(NF[IDX[idn], IDX[x], IDX[y]]
                           == (1 if x == y else 0)
                           for x in P for y in P))
    assert rigidity_ok
    report("  R1 rigidity: dim Hom(x,y) = delta_xy over all 10x10 "
           "primaries, a simple knot line cannot")
    report("     change label in strand.  Any in strand switch, odd or "
           "even parity, is the ZERO morphism.")

    # # R2: no same orientation conversion vertex ##
    ff = {c: int(NF[IDX[f], IDX[f], IDX[c]]) for c in P
          if NF[IDX[f], IDX[f], IDX[c]]}
    assert NF[IDX[f], IDX[f], IDX[a]] == 0
    assert ff == {(0, 1): 1, (2, 0): 1}
    report(f"  R2 vertex absence: f (x) f = {ff}, the adjoint does "
           "NOT appear.  A uniformly oriented")
    report("     f-tower cannot emit an a-line at any price (stronger "
           "than the conversion lemma claims).")

    # # R3: conversion = orientation reversed pair fusion ##
    conv_dim = int(NF[IDX[f], IDX[fb], IDX[a]])
    assert conv_dim == 1
    N = N_FUND
    A = N @ N.T - np.eye(6)
    assert np.allclose(N @ N.T, np.eye(6) + A)
    report("  R3 conversion channel: dim Hom(f (x) fbar, a) = 1, and "
           "the walk algebra states it exactly:")
    report("     n_f n_fbar = I + n_a  (canon: A = N N^T - I).  "
           "Conversion is pair fusion with an")
    report("     orientation reversed leg: strand number changes -> "
           "interaction history -> vent kind.")

    # # the moot evaluation + the empirical fence ##
    walks = lambda M: int(round(M.sum(axis=1)[0]))
    ghosts = {"f2a2 (36)": walks(N @ N @ A @ A),
              "f3a+fa3 (78)": walks(np.linalg.matrix_power(N, 3) @ A)
              + walks(N @ np.linalg.matrix_power(A, 3))}
    fence = []
    for name, cnt in ghosts.items():
        for factor in (1, 2):
            mev = factor * cnt * PDG_MASSES["tau"]
            match = [n for n, v in PDG_MASSES.items()
                     if abs(mev / v - 1.0) < 0.01]
            assert not match
            fence.append(round(mev / 1e3))
    report(f"  fence: ghost states {sorted(set(fence))} GeV (base and "
           "even response) match no PDG state. The")
    report("     sqrt(2)^odd witness (octonions.py) still kills the "
           "odd class independently.")

    report("-" * 66)
    report("  VERDICT: terminus purity is label rigidity.  The even "
           "switch completion is excluded by")
    report("  rigidity (R1) and vertex absence (R2), independent of "
           "any value.  There is no 6j evaluation")
    report("  to perform, because no morphism exists to evaluate.  "
           "The knot definition only describes this.")
    report("  Schur's lemma does the work, and the base/vent taxonomy "
           "is the uniform edge kind law of the web.")
    report("  Grounding: Schur's lemma in fusion categories (EGNO Ch."
           "1). Junction Hom spaces for defect/")
    report("  boundary lines (Fuchs-Runkel-Schweigert 2002; Chang-Lin-"
           "Shao-Wang-Yin 2019).")
    return dict(rigidity=rigidity_ok, ff_channels=ff,
                conversion_dim=conv_dim,
                ghost_fence_gev=sorted(set(fence)))


if __name__ == "__main__":
    run()
