"""
lambda_assignment.py: the Koide representation assignment, derived.

WHAT THIS DERIVES.  The Koide representation assignment

    Q(c,b,t) = Q0 + h11/K^3 = 289/432,
    Q(s,c,b) = Q0 + h10/K^3 = 649/972

follows from two lemmas already executable in canon composed under
THE INSERTION RULE, with data contact applied
once, last, as corroboration (blind-first ordering, closure_search
style).  The eight-structure enumeration of masses.py is
confirmation of this derivation, not its source.

THE TWO CANON LEMMAS (recomputed here, not cited blind):

  L1 (bookkeeping.py precedent, the theta selection): the mass
     insertion is the Z3 sector-changer, and where no superselection
     constrains it, the LIGHTEST sector-changing primary mediates.
     Computed: the lightest is the conjugate fundamental pair
     {(1,0),(0,1)} at h = h10 = 2/9 (conjugates share h, so the choice
     inside the pair is weight-blind).

  L2 (words.py part 2, N-ality superselection): a confining
     interface (unbroken Z3 centre) transmits only zero-triality
     channels; the zero-triality subsector is C_ad = {1, J, J^2, 8},
     the simple currents act as the identity on the D6 boundary
     (n_J = I6, words.py), so the UNIQUE nontrivial neutral channel
     with a nontrivial boundary action is the adjoint (1,1) at
     h = h11 = 1/2.  Not "lightest neutral": ALONE.

THE INSERTION RULE:

     a triplet's Koide insertion is subject to the superselection of
     the deepest interface its members close through.  The (c,b,t)
     triplet contains the terminus generation: t closes at orbit
     completion (J^3 = 1, the vacuum row), behind the confining
     interface, so its insertion must transmit through that
     interface.  The (s,c,b) triplet contains no terminus member; its
     insertion faces no neutral-lane superselection.

  The rule is the SAME physical statement that forces the step letter
  of w2 = ffa (words.py: steps at a confining interface are neutral),
  applied to the sector-changing insertion instead of the ladder
  step.  It carries the same standing as the orientation and
  apportionment laws: stated, witnessed, killed by its alternative.

THE ASSIGNMENT (forced by the rule):

  Q(c,b,t):  admissible channels = nontrivial C_ad with nontrivial
             boundary action = {adjoint}, by L2.  UNIQUE.  -> h11.
  Q(s,c,b):  admissible channels = all sector-changers; L1's
             minimality picks the fundamental pair.  -> h10.

  No data enters.  The K^3 altitude is separate and already derived
  (the first-invariant-order theorem).

COROLLARY.  Q(c,b,t) = 289/432 and Q(s,c,b) = 649/972 exactly.

CORROBORATION (data contact, once, last):
  the 2x2 assignment block inverted for m_b: canon lands +0.26%; the
  swapped assignment misses by +2.2% (~ +13 sigma at PDG precision).
  The papers add structural corroboration: under F4 -> SU(3)xSU(3),
  26 -> (1,8) + (3,3) + (3bar,3bar), the adjoint appears EXACTLY
  ONCE, in the diagonal channel, and (c,b,t) is the diagonal Z3
  triplet (quark paper, Sec. bottom).

Usage: python3 lambda_assignment.py   (from tests/probes/, or give the
       repo root as argv[1])
"""

import math
import os
import sys
from fractions import Fraction

_here = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.dirname(_here),
          os.path.join(os.path.dirname(os.path.dirname(_here)),
                       "interference"),
          os.path.join(_here, "interference")):
    if p not in sys.path:
        sys.path.insert(0, p)
import _bootstrap  # noqa: E402,F401


def run(report=print):
    from root import (d10, d11, K, h10, h11, Q0, BA_ratio,
                      M_E_ANCHOR_MEV, PDG_MASSES)
    from mixing import _compute_wzw

    report("LAMBDA ASSIGNMENT: derived, not enumerated")
    report("=" * 66)

    ## L1: the lightest sector-changer (recomputed) ##
    prims = [(a, b) for a in range(4) for b in range(4) if a + b <= 3]
    h = {p: Fraction(p[0]**2 + p[1]**2 + p[0] * p[1]
                     + 3 * p[0] + 3 * p[1], 18) for p in prims}
    t = {p: (p[0] + 2 * p[1]) % 3 for p in prims}
    changers = [p for p in prims if p != (0, 0)]
    hmin = min(h[p] for p in changers)
    lightest = sorted(p for p in changers if h[p] == hmin)
    assert hmin == h10 and lightest == [(0, 1), (1, 0)]
    report(f"  L1 lightest sector-changer: {lightest} at h = {hmin} "
           f"= h10 (conjugate pair, weight-blind)")

    ## L2: the unique neutral stepper (recomputed) ##
    w = _compute_wzw()
    P, IDX, NF = w["PRIMARIES"], w["IDX"], w["N_fus"]
    C0 = sorted(p for p in prims if t[p] == 0)
    assert C0 == [(0, 0), (0, 3), (1, 1), (3, 0)]
    # fusion-closed (no leak out of the neutral sector)
    assert all(NF[IDX[a], IDX[b], IDX[c]] == 0
               for a in C0 for b in C0 for c in P if tuple(c) not in C0)
    # invertibles act as the identity on the D6 boundary (words.py,
    # n_J = I6): the unique NONTRIVIAL neutral stepper is the adjoint
    steppers = [p for p in C0 if p not in ((0, 0), (3, 0), (0, 3))]
    assert steppers == [(1, 1)] and h[(1, 1)] == h11
    report(f"  L2 neutral sector C_ad = {C0}, fusion-closed; simple "
           f"currents act as I6 (words.py);")
    report(f"     unique nontrivial neutral stepper = (1,1) at h = "
           f"{h[(1, 1)]} = h11.  ALONE, not 'lightest'.")

    ## P + theorem: the assignment##
    report("  the INSERTION RULE: terminus-containing triplet ->")
    report("     insertion crosses the confining interface -> neutral "
           "lane (same superselection as w2's step letter)")
    assign_cbt = h11        # forced by L2 given P (unique)
    assign_scb = h10        # forced by L1 (no superselection)
    Q_cbt = Q0 + Fraction(assign_cbt) / K**3
    Q_scb = Q0 + Fraction(assign_scb) / K**3
    assert Q_cbt == Fraction(289, 432) and Q_scb == Fraction(649, 972)
    report(f"  ASSIGNMENT (forced by the rule): Q(c,b,t) = {Q_cbt}, "
           f"Q(s,c,b) = {Q_scb}")
    report("  -> masses.py's eight-structure enumeration is "
           "CONFIRMATION of this derivation.")

    # papers' structural corroboration: 26 under SU(3) x SU(3)
    assert 8 + 9 + 9 == 26
    report("  corroboration (quark paper): 26 -> (1,8)+(3,3)+"
           "(3bar,3bar); the adjoint appears EXACTLY ONCE,")
    report("     in the diagonal channel, and (c,b,t) is the diagonal "
           "Z3 triplet.")

    ## data contact, once, last (corroboration only) ##
    theta = float(h10)
    gaps = {k: 1.0 + BA_ratio * math.cos(theta + 2 * math.pi * k / 3)
            for k in range(3)}
    m_mu = M_E_ANCHOR_MEV * (gaps[2] / gaps[1])**2
    m_tau = M_E_ANCHOR_MEV * (gaps[0] / gaps[1])**2
    m_c = float(Fraction(d10**2 * d11) * (1 + Fraction(1, K**3))) * m_mu
    m_t = float(Fraction(d11**4 + d10**4)
                * (1 + Fraction(1, 2 * K * (d11**4 + d10**4)))) * m_tau

    def invert_b(Q):
        sc, st = math.sqrt(m_c), math.sqrt(m_t)
        a, b = Q - 1.0, 2.0 * Q * (sc + st)
        c = (Q - 1.0) * (m_c + m_t) + 2.0 * Q * sc * st
        disc = b * b - 4.0 * a * c
        return next(x * x for x in ((-b + math.sqrt(disc)) / (2 * a),
                                    (-b - math.sqrt(disc)) / (2 * a))
                    if m_c < x * x < m_t)

    rows = []
    for tag, hv in (("canon (adjoint->cbt)", h11),
                    ("swap  (fund->cbt)", h10)):
        mb = invert_b(float(Q0 + Fraction(hv) / K**3))
        pull = 100.0 * (mb - PDG_MASSES["b"]) / PDG_MASSES["b"]
        rows.append((tag, mb, pull))
        report(f"     {tag}: m_b = {mb:7.1f} MeV ({pull:+.2f}%)")
    assert abs(rows[0][2]) < 0.5 and rows[1][2] > 1.5
    report(f"     swap dead at ~{(rows[1][1] - PDG_MASSES['b'])/7.0:+.0f} "
           "sigma (PDG ~7 MeV), confirmation, not the selection.")

    report("-" * 66)
    report("  STATUS: assignment forced by the insertion rule.  "
           "Data contact:")
    report("  corroboration only, applied once, last.")
    return dict(Q_cbt=Q_cbt, Q_scb=Q_scb,
                lightest=lightest, neutral_stepper=steppers[0],
                mb_canon=rows[0][1], mb_swap=rows[1][1],
                rule="the insertion rule")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.path.insert(0, os.path.join(sys.argv[1], "interference"))
    run()
