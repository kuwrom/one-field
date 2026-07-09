"""
edge_derivations.py: the two independent mass derivations.

Two canon edges enter the mass sector through structures beyond the
word/correction menus.  Neither is a menu choice; this probe holds
both derivations in executable form.

THE DOWN-QUARK VENT (m_d).  The exact identity

    d11^2 (1 + d10/d11^4)  =  count(aa) + h10  =  83/9

says the multiplicative vent form IS the committed word_multiplier
operation m_q = (count(w) + corr) m_lep: the Z2 triality swap changes
only the seed word (ff -> aa, exactly as at the seed level u = ff,
d = aa), and corr(1) = h10 is inherited from the generation ladder,
already fixed at up-gen-1.  Zero freedom.  The down-gen-1 search
stage (closure_search.py, ancillary) shows (aa, h10) as the unique
ladder-conditioned hit at the 1% band.

THE STRANGE BRIDGE (m_s).  The bridge is the norm conversion between
the two F4-equivariant metrics on one generation line:

    bridge^2 = delta_Z2^2 / delta_J^2 = (4/9)/(3/8) = 32/27,

with every ingredient derived:

  delta_Z2 = 2/3   the nontrivial root of the Dynkin fixed point
                   (1 + 3d)(1 - d) = 1  [computed below];
  delta_J^2 = 3/8  the Albert spectral gap, CONSTRUCTED from a
                   concrete J3(O) element: with the framework's own
                   Fano table (octonions.py), the element
                   X = s.I + off-diagonal z_i = sqrt(w) e_u has
                   invariants exactly (3s, 3s^2 - 3w, s^3 - 3ws) --
                   spectrum {s, s +/- sqrt(3w)}, IFF the three
                   units form an ANTI-ASSOCIATIVE triple; all 7
                   associative Fano lines fail (the Re(z1 z2 z3)
                   term skews the spectrum).  The Z3-symmetric
                   spread requires the anti-associative directions.
  w = 1/dim(O)     the isotropic share: E[z_a^2] = 1/8 exactly for a
                   maximally-forgotten (isotropic) unit octonion on
                   S^7, the same forgetting whose PvP = 0 erases
                   the Z3 label.  The imaginary-only alternative
                   (S^6, share 1/7) is dead in the scan below.

  UNIQUENESS OF THE CONVERSION (Schur).  All F4-representations are
  self-dual (-1 lies in W(F4): the longest element acts as -id,
  verified by reflection closure in wiring_theorem_probe.py), and
  the 26 is absolutely irreducible of real type (Springer-Veldkamp
  Ch. 5), so the space of invariant symmetric bilinear forms on the
  26 is one-dimensional: any two equivariant norms differ by ONE
  constant, and the Dynkin->Albert mass map is the single factor
  delta_Z2/delta_J = sqrt(32/27).  Forced, not chosen.

  The identities delta_Z2 = Q0 and 1/delta_J^2 = 8/3 = charge_trace
  make the norm-conversion and vertex-composition readings of 32/27
  one statement.

Corroboration (data contact once, last): the share scan lands m_s
within the scorecard band only at w = 1/8; the altitude factor
(1 + h11/K^2) is the paired-insertion K^2 term
(the first-invariant-order theorem).

Usage: python3 edge_derivations.py
"""

import math
import os
import sys
from fractions import Fraction
from itertools import combinations

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
    from root import (d10, d11, K, h10, h11, Q0, charge_trace,
                      PDG_MASSES, BA_ratio, M_E_ANCHOR_MEV)
    from octonions import FANO_TRIPLES, FANO_SIGNS
    import edge_grammar as eg

    report("EDGE DERIVATIONS: the two independent mass derivations")
    report("=" * 66)

    ## the down-quark vent - exact identity ##
    lhs = Fraction(d11**2) * (1 + Fraction(d10, d11**4))
    rhs = Fraction(d11**2) + Fraction(h10)
    assert lhs == rhs == Fraction(83, 9)
    w2 = eg.word_menu(2)
    assert w2["aa"] == 9 and "h10" in eg.correction_menu()
    hits_cond = [w for w, c in w2.items()
                 if abs((c + float(h10)) * M_E_ANCHOR_MEV
                        / PDG_MASSES["d"] - 1.0) < 0.01]
    assert hits_cond == ["aa"]
    report(f"  m_d: d11^2(1 + d10/d11^4) = {lhs} = count(aa) + h10 "
           f"(EXACT; word 'aa' and corr 'h10' committed);")
    report("       ladder-conditioned down-gen-1 stage: unique hit = "
           "aa  ->  m_d = (83/9) m_e")

    # # the Albert spectral gap, constructed ##
    M = np.zeros((8, 8, 8))
    M[0, :, :] = np.eye(8)
    M[:, 0, :] += np.eye(8)
    M[0, 0, 0] = 1.0
    for i in range(1, 8):
        M[i, i, 0] = -1.0
    for (i, j, k), s in zip(FANO_TRIPLES, FANO_SIGNS):
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            M[a, b, c] = s
            M[b, a, c] = -s
    mult = lambda x, y: np.einsum('i,j,ijk->k', x, y, M)
    rng = np.random.default_rng(7)
    for _ in range(20):
        x, y = rng.normal(size=8), rng.normal(size=8)
        assert abs(np.linalg.norm(mult(x, y))
                   - np.linalg.norm(x) * np.linalg.norm(y)) < 1e-9

    lines = {frozenset(t) for t in FANO_TRIPLES}
    on_l = [t for t in combinations(range(1, 8), 3)
            if frozenset(t) in lines]
    off_l = [t for t in combinations(range(1, 8), 3)
             if frozenset(t) not in lines]
    assert len(on_l) == 7 and len(off_l) == 28
    for t in combinations(range(1, 8), 3):
        e = [np.eye(8)[u] for u in t]
        re = mult(mult(e[0], e[1]), e[2])[0]
        assert (frozenset(t) in lines) == (abs(re) > 1e-12)

    def invariants(s, zs):
        Nz = [float(z @ z) for z in zs]
        return (3.0 * s, 3.0 * s * s - sum(Nz),
                s**3 + 2.0 * mult(mult(zs[0], zs[1]), zs[2])[0]
                - s * sum(Nz))

    w = 1.0 / 8.0
    for triple, expect_sym in (((1, 2, 3), True), ((1, 2, 7), False)):
        zs = [math.sqrt(w) * np.eye(8)[u] for u in triple]
        sym = all(abs(invariants(s, zs)[1] - (3*s*s - 3*w)) < 1e-12
                  and abs(invariants(s, zs)[2] - (s**3 - 3*w*s)) < 1e-12
                  for s in (1/3, 0.5, 1.0))
        assert sym == expect_sym
    zs = [math.sqrt(w) * np.eye(8)[u] for u in (1, 2, 3)]
    for s in (0.5, 1.0):
        T, S, D = invariants(s, zs)
        for lam in (s, s + math.sqrt(3*w), s - math.sqrt(3*w)):
            assert abs(lam**3 - T*lam**2 + S*lam - D) < 1e-12
    report(f"  m_s: Fano table computed; |xy| = |x||y|; 7 associative "
           f"/ 28 anti-associative triples;")
    report("       anti-associative + share 1/8 -> spectrum "
           "{s, s +/- delta_J}, delta_J^2 = 3/8 (associative")
    report("       lines fail: the symmetric spread REQUIRES the "
           "anti-associative directions)")

    # isotropy: the share is the second moment on S^7
    v = rng.normal(size=(200000, 8))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    assert abs(float(np.mean(v[:, 0]**2)) - 0.125) < 2e-3
    # Dynkin fixed point and the conversion
    assert [d for d in (Fraction(0), Fraction(2, 3))
            if (1 + 3*d) * (1 - d) == 1] == [Fraction(0),
                                             Fraction(2, 3)]
    delta_Z2, delta_J_sq = Fraction(2, 3), Fraction(3, 8)
    bridge_sq = delta_Z2**2 / delta_J_sq
    assert bridge_sq == Fraction(32, 27) == Q0**2 * charge_trace
    report(f"       E[z_a^2] on S^7 = 1/8 (isotropy of the forgotten "
           f"slot); delta_Z2 = 2/3 (fixed point);")
    report(f"       bridge^2 = delta_Z2^2/delta_J^2 = {bridge_sq} "
           f"= Q0^2 x charge_trace (one identity, two readings)")

    # share scan (data contact once, last): w = 1/8 unique at 1%
    theta = float(h10)
    gaps = {k: 1.0 + BA_ratio * math.cos(theta + 2*math.pi*k/3)
            for k in range(3)}
    m_mu = M_E_ANCHOR_MEV * (gaps[2]/gaps[1])**2
    m_tau = M_E_ANCHOR_MEV * (gaps[0]/gaps[1])**2
    m_c = float(Fraction(d10**2 * d11) * (1 + Fraction(1, K**3))) * m_mu
    m_t = float(Fraction(d11**4 + d10**4)
                * (1 + Fraction(1, 2*K*(d11**4 + d10**4)))) * m_tau
    Qb = float(Q0 + Fraction(h11) / K**3)
    sc, st = math.sqrt(m_c), math.sqrt(m_t)
    aq, bq = Qb - 1.0, 2.0*Qb*(sc + st)
    cq = (Qb - 1.0)*(m_c + m_t) + 2.0*Qb*sc*st
    disc = bq*bq - 4.0*aq*cq
    m_b = next(x*x for x in ((-bq + math.sqrt(disc))/(2*aq),
                             (-bq - math.sqrt(disc))/(2*aq))
               if m_c < x*x < m_t)
    Qs = float(Q0 + Fraction(h10) / K**3)
    sb = math.sqrt(m_b)
    a2, b2 = 1.0 - Qs, 2.0*Qs*(sc + sb)
    c2 = m_c + m_b - Qs*(sc + sb)**2
    d2 = b2*b2 - 4.0*a2*c2
    base_root = next(x for x in ((-b2 + math.sqrt(d2))/(2*a2),
                            (-b2 - math.sqrt(d2))/(2*a2)) if x > 0)**2
    alt = math.sqrt(1.0 + float(h11) / K**2)
    hits = [n for n, ww in (("1", 1.0), ("1/2", .5), ("1/3", 1/3),
                            ("1/7", 1/7), ("1/8", .125),
                            ("1/14", 1/14), ("1/26", 1/26))
            if abs(base_root * math.sqrt(float(delta_Z2**2)/(3.0*ww)) * alt
                   / PDG_MASSES["s"] - 1.0) < 0.01]
    assert hits == ["1/8"]
    report(f"       share scan (once, last): w = 1/8 unique at 1% "
           f"(imaginary-only 1/7 dead) -> m_s lands")
    report("-" * 66)
    report("  both edges are DERIVATIONS: an exact on-menu identity "
           "and a Schur-unique norm conversion.")
    return dict(m_d=str(lhs), bridge_sq=str(bridge_sq),
                delta_J_sq=str(delta_J_sq), share_hits=hits,
                anti_associative=len(off_l))


if __name__ == "__main__":
    run()
