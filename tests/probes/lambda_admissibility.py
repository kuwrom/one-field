"""
lambda_admissibility.py: the Cabibbo reading, from its declared menu to one entry.

THE MENU. mixing.py reads lambda as the Sugawara weight h10 carrying one
echo, lambda = h10 (1 + x / K^n), with the inserted weight x in
{h10, h11, Q0, 1} and the altitude K^n in {K, K^2, K^3}: twelve entries.
This probe cuts the menu down with lemmas the engine already holds and one
rule stated in the Coherence Gate paper. It uses no measured value.

STATEMENTS (computed below).

  A1 (order).  The insertion is the generation shift S of the circulant.
     tr S = tr S^2 = 0 and tr S^3 = 3, so a single insertion first returns
     to the identity at order three, altitude K^3.  A paired insertion
     S S^dagger is the identity at order two, altitude K^2.  No insertion
     enters at K.
  A2 (weights).  The inserted x is the conformal weight of the channel the
     insertion closes through, a weight of SU(3)_3: {0, 2/9, 1/2, 5/9, 8/9,
     1}.  Q0 = 2/3 is the Koide constant, not a weight.
  A3 (the unit).  x = 1 is the weight of the simple current J = (3,0).  Its
     boundary action on the D6 boundary is the identity, n_J = 1, built by
     fusion from the fundamental nimrep, so J closes no transition between
     generations.
  A4 (the paired channel).  A paired insertion closes through
     f x fbar = 1 + 8, so its nontrivial channel is the adjoint, x = h11.
  A5 (the overlap).  The mixing matrix is the overlap V = U_u^dagger U_d of
     two rotations, one entering conjugated.  By the Hermitian pair rule of
     the Octavian Coherence Gate paper, Sec. X ("the off diagonal coupling
     appears once for each of the two orientations S, S^dagger, giving
     B = 2|c_1|"),
     its insertion is the pair S S^dagger.  Single insertions do not give a
     reading of V.

RESULT. A1 to A4 leave three entries: paired h11 / K^2, single h10 / K^3 and
single h11 / K^3.  A5 leaves one, lambda = h10 (1 + h11/K^2) = (2/9)(73/72)
= 73/324. This is the same paired insertion factor 73/72 that root.py
applies to bridge^2.  The data check it (tests/mixing/test_mixing_probes.py).

Usage: python3 tests/probes/lambda_admissibility.py     (instant)
"""

from fractions import Fraction

import numpy as np

from _bootstrap import *  # noqa: F401,F403
from interference import mixing
from interference.root import K, Q0, d10, d11, h10, h11


def weights():
    """Conformal weights h = C2 / (k + 3) of the ten SU(3)_3 primaries."""
    out = {}
    for a, b in mixing._compute_wzw()["PRIMARIES"]:
        c2 = Fraction(a * a + b * b + a * b + 3 * a + 3 * b, 3)
        out[(a, b)] = c2 / (3 + 3)
    return out


def order_lemma():
    """First return to the identity: single S at order 3, paired S S^dagger at order 2."""
    S = np.roll(np.eye(3), 1, axis=0)
    traces = [round(float(np.trace(np.linalg.matrix_power(S, n)))) for n in (1, 2, 3)]
    paired_identity = bool(np.allclose(S @ S.T, np.eye(3)))
    single = min(n for n in (1, 2, 3) if traces[n - 1] != 0)
    return dict(traces=traces, single=single, paired=2 if paired_identity else None)


def current_action():
    """n_J for J = (3,0) on the D6 boundary, by fusion from the fundamental nimrep."""
    d6 = mixing._compute_d6_boundary()
    n10 = np.asarray(d6["n10"], dtype=float)
    n01, eye = n10.T, np.eye(len(n10))
    n11 = n10 @ n01 - eye                 # (1,0)(0,1) = (0,0) + (1,1)
    n20 = n10 @ n10 - n01                 # (1,0)(1,0) = (2,0) + (0,1)
    n30 = n20 @ n10 - n11                 # (2,0)(1,0) = (3,0) + (1,1)
    return bool(np.allclose(n30, eye)), bool(np.allclose(n11, np.asarray(d6["n11"], dtype=float)))


def paired_channel():
    """The nontrivial channel of f x fbar."""
    w = mixing._compute_wzw()
    I, N = w["IDX"], w["N_fus"]
    f, fb = (1, 0), (0, 1)
    out = [p for p in w["PRIMARIES"] if round(float(N[I[f]][I[fb], I[p]])) > 0]
    return sorted(out)


def menu():
    xs = {"h10": Fraction(h10), "h11": Fraction(h11), "Q0": Fraction(Q0), "1": Fraction(1)}
    return [(name, x, n) for name, x in xs.items() for n in (1, 2, 3)]


def run(report=print):
    entries = menu()
    w = weights()
    order = order_lemma()
    j_identity, n11_consistent = current_action()
    channel = paired_channel()
    assert order["traces"] == [0, 0, 3] and order["single"] == 3 and order["paired"] == 2
    assert j_identity and n11_consistent
    assert channel == [(0, 0), (1, 1)] and w[(1, 1)] == Fraction(h11) and w[(1, 0)] == Fraction(h10)
    weight_set = set(w.values())

    a1 = [e for e in entries if e[2] in (order["single"], order["paired"])]
    a2 = [e for e in a1 if e[1] in weight_set]
    a3 = [e for e in a2 if e[1] != w[(3, 0)]]
    adj = w[(1, 1)]
    a4 = [(nm, x, n, "paired" if n == order["paired"] else "single") for nm, x, n in a3
          if not (n == order["paired"] and x != adj)]
    a5 = [e for e in a4 if e[3] == "paired"]
    assert len(entries) == 12 and len(a4) == 3 and len(a5) == 1
    name, x, n, kind = a5[0]
    lam = Fraction(h10) * (1 + x / K**n)
    assert lam == Fraction(73, 324) and 1 + x / K**n == Fraction(73, 72)

    report(f"A1 order: tr S^n = {order['traces']}. Single insertions at K^{order['single']}, "
           f"paired at K^{order['paired']}: {len(entries)} -> {len(a1)} entries")
    report(f"A2 weights {sorted(weight_set)}; Q0 = {Fraction(Q0)} is not one: -> {len(a2)}")
    report(f"A3 n_J is the identity on the D6 boundary: x = 1 closes no transition: -> {len(a3)}")
    report(f"A4 f x fbar = {channel}: a paired insertion carries h11 = {adj}: -> {len(a4)}: "
           + ", ".join(f"{k} {nm}/K^{m}" for nm, _, m, k in a4))
    report(f"A5 the overlap V = U_u^dagger U_d is a Hermitian pair: -> {len(a5)}: "
           f"lambda = h10 (1 + {name}/K^{n}) = {lam}")
    return dict(menu=len(entries), after=[len(a1), len(a2), len(a3), len(a4), len(a5)],
                admissible=[(nm, m, k) for nm, _, m, k in a4], survivor=(name, n), lam=lam)


if __name__ == "__main__":
    run()
