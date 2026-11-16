"""The PMNS sign as two framework bits, both fixed by data blind gates.

THE QUESTION. The sign of the matter excess only decides which kind is
called matter, because one orientation bit flips every CP-odd sign
together.  The physical content is the
RELATIVE sign of the lepton and quark CP phases.  This probe computes
which framework readings fix it without data.

THE TWO BITS. The Standard Model writes V_CKM = U_u^dag U_d and
U_PMNS = U_e^dag U_nu.  The engine builds V_CKM directly from its
four readings (mixing.py), with e^{-i delta_q} at (1,3) for the
orientation that gives eta bar = +pi/9 = arg S[(1,0),(1,0)].  The
lepton phase sits in U_e (U_TBM is real).
  q  which quark rotation carries the phase: U_d carries e^{-i delta}
     at (1,3); U_u, entering daggered, carries e^{+i delta}.  CKM
     observables read only V, so q is read from the charges (P4).
  p  which pairing the lepton sector reads its orientation through:
     the diagonal pairing (fundamental, fundamental), arg S = +pi/9, or
     its own conjugation invariant Z_C (fundamental, antifundamental),
     arg S[(1,0),(0,1)] = -pi/9 = conj, since S[a, b*] = conj S[a, b].
Take the phase carrying rotation of each sector to carry the same
orientation in the same slot, up to p's conjugation. Then U_e carries
e^{s_q s_p i a} at (1,3), a = arccos(-lambda). This matching rule is
a reading, stated here.

P1. S[a, b*] = conj S[a, b] for every primary of SU(3)_3.
P2. The invariant delta for the four (q, p) rows: 275.74 deg for
    (U_d, diagonal) and (U_u, Z_C); 84.26 deg for the other two.  Only
    the product of the two bits is physical.
P3. The lepton sector reads Z_C for its angles (the conjugation orbit
    TBM block, mixing.py). One sector reads one invariant, so this
    forces p = Z_C.  gates() builds mixing.py's generation block under each
    pairing: Z_C gives d11 carriers with a democratic adjoint block,
    the diagonal pairing gives no paired orbit and no block.  With
    p = Z_C the two remaining rows split by q: 275.74 deg for U_u,
    84.26 deg for U_d.
P4. The charge dictionary (README, Quantum numbers; root.py) reads
    Q_u = d10/d11 with d10 = S[0,(1,0)]/S[0,0] = 2, the quantum
    dimension of the fundamental primary, and Q_d = 1/d11, which reads
    no (1,0) datum.  The fundamental primary carries the up-type
    sector.  eta bar is arg S[(1,0),(1,0)], the same primary's
    oriented phase, so the up-type rotation carries the CP phase:
    q = U_u. This is the charge dictionary reading (IH). The dictionary
    is made for the charges. gates() takes the phase carrying rotation
    from it: the only sector whose charge depends on (1,0).

GATES. gates() computes P3 and P4 from the modular data alone: of the
four rows exactly one passes both, (U_u, Z_C), delta = 275.74 deg.
The NuFIT 3-sigma range checks it (sign_exhaustion_theorem.py reaches
the same row from the data).

RESULT. P1 gives the conjugation of the pairing phase.  P3 fixes
p = Z_C from the lepton sector's own invariant.  P4 fixes q = U_u
from the charge dictionary (IH). The matching rule places the phase.
The four rows of P2 reduce to one, delta = 275.74 deg, and the NuFIT
range is the check. Every CKM observable depends on V alone. So the
relative CP sign is the only observable that depends on q, just as G
is the only observable that depends on the boundary reading. P4 has
the same form as P3: one primary, one sector.

Usage: python3 tests/probes/pmns_sign_bits.py     (instant)
"""

import cmath
from fractions import Fraction
import math

import numpy as np

from _bootstrap import *  # noqa: F401,F403
from interference import mixing
from interference.root import BA_ratio, alpha_phys, d10, d11, h10

TBM = np.array([[math.sqrt(2 / 3), 1 / math.sqrt(3), 0],
                [-1 / math.sqrt(6), 1 / math.sqrt(3), 1 / math.sqrt(2)],
                [1 / math.sqrt(6), -1 / math.sqrt(3), 1 / math.sqrt(2)]])


def _R13(t, d):
    c, s, e = math.cos(t), math.sin(t), cmath.exp(1j * d)
    return np.array([[c, 0, s * e], [0, 1, 0], [-s * e.conjugate(), 0, c]])


def _R23(t):
    c, s = math.cos(t), math.sin(t)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def invariant_delta(U):
    s13 = abs(U[0, 2])
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = abs(U[0, 1]) / c13, abs(U[1, 2]) / c13
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    J = float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))
    sd = J / (s12 * c12 * s23 * c23 * s13 * c13**2)
    cd = (abs(U[1, 0])**2 - s12**2 * c23**2 - c12**2 * s23**2 * s13**2) / (2 * s12 * c12 * s23 * c23 * s13)
    return math.degrees(math.atan2(sd, cd)) % 360


def rows():
    lam = float(h10) * (1.0 + 0.5 / (d10 * d11)**2)
    a = math.acos(-lam)

    def mf(k):
        gap = 1.0 + BA_ratio * math.cos(float(h10) + 2 * math.pi * k / 3)
        return gap**2 * (1.0 - (d11**3 / d10) * (alpha_phys / (2 * math.pi))**2 * gap)
    phi = math.sqrt(mf(1) / mf(2)) / math.sqrt(d10)
    out = {}
    for q, s_q in (("U_d", -1), ("U_u", +1)):
        for p, s_p in (("diagonal", +1), ("Z_C", -1)):
            U = np.conj((_R23(-phi) @ _R13(float(h10), s_q * s_p * a)).T) @ TBM
            out[(q, p)] = invariant_delta(U)
    return out


def identity():
    w = mixing._compute_wzw()
    S, I, P = w["S_mat"], w["IDX"], w["PRIMARIES"]
    worst = max(abs(S[I[x], I[(y[1], y[0])]] - np.conj(S[I[x], I[y]])) for x in P for y in P)
    f, fb = (1, 0), (0, 1)
    return worst, cmath.phase(S[I[f], I[f]]) / math.pi, cmath.phase(S[I[f], I[fb]]) / math.pi


def dictionary():
    w = mixing._compute_wzw()
    S, I = w["S_mat"], w["IDX"]
    d_fund = float(np.real(S[0, I[(1, 0)]] / S[0, 0]))
    return d_fund, d_fund / float(d11), 1.0 / float(d11)


SIGNS = {"U_d": -1, "U_u": +1, "diagonal": +1, "Z_C": -1}


def gates():
    """P3 and P4 as data blind gates on the four (q, p) rows.

    P3: the pairing the lepton sector reads its angles through.  The
    generation block of mixing.py is built from the orbits of a pairing.
    Under Z_C, (a, a*), three conjugate pairs carry a democratic adjoint
    block, the TBM block.  Under the diagonal pairing, (a, a), no orbit
    is a pair and no generation block forms.  Only Z_C gives d11 carriers.
    P4: the rotation that carries the CP phase.  eta bar is the phase of
    S[(1,0),(1,0)], so the phase rides the sector whose charge reads the
    primary (1,0): Q_u = d(1,0)/d(1,1) reads it, Q_d = 1/d(1,1) does not.
    """
    w = mixing._compute_wzw()
    S, I, P, N = w["S_mat"], w["IDX"], w["PRIMARIES"], w["N_fus"]
    qdim = {x: float(np.real(S[0, I[x]] / S[0, 0])) for x in P}

    def generation_block(pair_map):
        """mixing.py's construction under a pairing: orbits of the map, the
        noninvertible paired orbits as carriers, the adjoint's folded block."""
        orbits, seen = [], set()
        for x in P:
            if x not in seen:
                y = pair_map[x]
                orbits.append((x,) if y == x else (x, y))
                seen.update([x, y])
        Q = np.zeros((len(P), len(orbits)))
        for k, o in enumerate(orbits):
            for x in o:
                Q[I[x], k] = 1 / math.sqrt(len(o))
        gen = [k for k, o in enumerate(orbits) if len(o) == 2 and qdim[o[0]] > 1 + 1e-9]
        blk = (Q.T @ N[I[(1, 1)]] @ Q)[np.ix_(gen, gen)]
        return len(gen), bool(len(gen) and np.allclose(blk, blk[0, 0]))

    maps = {"Z_C": {x: (x[1], x[0]) for x in P}, "diagonal": {x: x for x in P}}
    blocks = {name: generation_block(m) for name, m in maps.items()}
    passing = [name for name, (n, dem) in blocks.items() if n == d11 and dem]
    assert len(passing) == 1
    pairing = passing[0]
    f, adj = (1, 0), (1, 1)
    assert abs(cmath.phase(S[I[f], I[f]]) - math.pi / 9) < 1e-12
    # Which primaries each charge reads, from root.py's forms: a charge
    # reads a primary when its form depends on that primary's dimension.
    import sympy as sp
    from interference import root
    d_f, d_adj = sp.symbols("d_f d_adj", positive=True)
    symbol_of = {d_f: f, d_adj: adj}
    reads = {rot: tuple(symbol_of[x] for x in sorted(sp.sympify(root.CHARGE_FORMS[q](d_f, d_adj)).free_symbols,
                                                      key=str))
             for rot, q in (("U_u", "Q_u"), ("U_d", "Q_d"))}
    assert root.Q_u == Fraction(2, 3) and root.Q_d == Fraction(1, 3)
    for q in ("Q_u", "Q_d"):                         # the S-matrix dimensions give the same charges
        assert abs(root.CHARGE_FORMS[q](qdim[f], qdim[adj]) - float(getattr(root, q))) < 1e-12
    rotation = [q for q, prim in reads.items() if f in prim]
    assert len(rotation) == 1
    survivors = [(q, p) for q in ("U_d", "U_u") for p in ("diagonal", "Z_C")
                 if p == pairing and q == rotation[0]]
    return dict(pairing=pairing, blocks=blocks, rotation=rotation[0], survivors=survivors,
                delta_sign=SIGNS[survivors[0][0]] * SIGNS[survivors[0][1]])


def run(report=print):
    worst, diag, zc = identity()
    r = rows()
    assert worst < 1e-12 and abs(diag - 1 / 9) < 1e-12 and abs(zc + 1 / 9) < 1e-12
    assert abs(r[("U_d", "diagonal")] - r[("U_u", "Z_C")]) < 1e-9
    assert abs(r[("U_u", "diagonal")] - r[("U_d", "Z_C")]) < 1e-9
    assert abs(r[("U_u", "Z_C")] - 275.7354) < 1e-3 and abs(r[("U_d", "Z_C")] - 84.2646) < 1e-3
    report(f"P1 S[a, b*] = conj S[a, b] to {worst:.1e}. Arg S[f, f] = {diag:+.4f} pi, "
           f"arg S[f, fbar] = {zc:+.4f} pi.")
    report("P2 invariant delta: " + "; ".join(f"{q}, {p}: {d:.2f} deg" for (q, p), d in r.items())
           + ".  Only the product of the two bits is physical.")
    report(f"P3 the lepton sector reads Z_C for its angles, so p = Z_C.  With p = Z_C the rows "
           f"split by q: U_u gives {r[('U_u', 'Z_C')]:.2f} deg, U_d gives {r[('U_d', 'Z_C')]:.2f} deg.")
    d_fund, Qu, Qd = dictionary()
    assert abs(d_fund - float(d10)) < 1e-12 and abs(Qu - 2 / 3) < 1e-12 and abs(Qd - 1 / 3) < 1e-12
    report(f"P4 d(1,0) = {d_fund:.6f} = d10; Q_u = d(1,0)/d11 = {Qu:.6f} reads the fundamental, "
           f"Q_d = 1/d11 = {Qd:.6f} does not: the fundamental carries the up-type sector and its "
           f"oriented phase is eta bar, so q = U_u and delta = {r[('U_u', 'Z_C')]:.2f} deg.")
    g = gates()
    assert g["survivors"] == [("U_u", "Z_C")]
    report(f"GATES P3 and P4, no data: the lepton pairing is {g['pairing']}, the phase sits on "
           f"{g['rotation']}. Of four rows one survives, delta = {r[g['survivors'][0]]:.2f} deg.")
    return dict(identity=worst, rows={f"{q}|{p}": d for (q, p), d in r.items()},
                d_fund=d_fund, Q_u=Qu, Q_d=Qd,
                survivors=[f"{q}|{p}" for q, p in g["survivors"]], delta_sign=g["delta_sign"])


if __name__ == "__main__":
    run()
