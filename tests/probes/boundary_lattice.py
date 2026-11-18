"""The boundary reading c on the whole lattice of shares: the one reading that uses no label.

THE QUESTION. boundary_whole.py fixes c = 1 on a menu of twelve: the
pointed orbit or the whole, the local response, its remainder or the
whole, one relative channel or both. The gates as wired allow c between
0.656 and 4.11 with no measured G (tests/audit/g_path_trials.py G10),
and the menu does not fill that window. Does the engine's coupling
definition leave a continuum for c, and what selects one reading on the
full lattice with no measured G?

PREMISES. P1 the engine's channel coupling: "closure probability times
the normalized dimension one current response" (interference/wiring.py),
with the second moment of the harmonic where the common mode enters
(boundary_whole.py P2). P2 each factor is a partial sum of a list: the
squares of the S-matrix's first row over the ten primaries, the two
parts of the unit response, the powers of the two relative channels.
P3 the horizon: the common mode past it carries no names
(interference/__init__.py, THE HORIZON). P4 the gates S and M on the
one list (interference/boundary.py), which use no measured G.

STATEMENTS (the numbers below check each one).

  L1 (the lattice, the reading's zero edge term). A reading is a
     nonempty subset of each list. The 1023 subsets of the primaries
     give 36 probabilities, k/36. With the three response readings and
     the two moment readings, the coefficient G M_b^2 has 216 readings
     and 162 values. Each factor is at most one, so every zero edge
     reading has c >= 1, and c = 1 is attained once. The window the
     gates leave below one contains no zero edge reading. The carrier
     edge of L7 adds to the reading and is not an entry of the lists.
  L2 (no continuum). The reading next to c = 1 is 36/35, with one
     simple current left out. It gives G lower by 2.8 percent, 1236
     standard deviations from the recommended value. So the zero edge
     term is discrete.
  L3 (the gates do not select). 71 of the 162 values pass S and M as
     wired at every contact, and four pass S under all six of its
     conventions: c = 1, 36/35, 18/17 and 12/11.
  L4 (symmetry does not select). The permutations of the ten primaries
     that keep the S-matrix and the conformal weights are two, the
     identity and charge conjugation. Conjugation keeps 22 values of
     c, and all four that pass S under every convention are among
     them: 1, 36/35, 18/17 and 12/11.
  L5 (the closures the common mode reads are closed under fusion).
     Closures compose by fusion, so a set of closures that the common
     mode reads contains the vacuum and every fusion product of its
     members. There are four fusion closed sets: the vacuum alone, the
     pointed orbit, the triality zero primaries and the whole, with
     c = 36, 12, 3 and 1. The whole without the pointed orbit contains
     no vacuum and is none of them. Every closure reads one common
     mode (root.py, THE COMMON MODE). The
     fundamental, the lightest triality changer, which sets theta
     (masses.py Step 3), is a closure. Of the four sets, only the whole
     contains it. FAILS for the other three: c = 36 and 12 fail S and
     M, and c = 3 passes them as wired but contains no triality
     changer.
  L5b (the input is step 3). The engine's closures sit in channels by
     their weights. The leptons' phase theta = h10 and the (s, c, b)
     backreaction h10/K^3 belong to the fundamental, and the (c, b, t)
     backreaction h11/K^3 belongs to the adjoint. The triality zero set
     contains the adjoint and not the fundamental. A common mode that
     read only that set would read the (c, b, t) triplet and neither
     the leptons nor (s, c, b). A closure that the common mode does not
     read would not move by any fraction of a bias. The ratio of its
     clock to the clock of a read closure would record this, which is
     the failure gate U names (gates_as_records.py U). The statement
     that every closure reads the one common mode is step 3 of the
     step table, and the exclusion of c = 3 rests on it. Of the four
     fusion closed sets, only the whole contains every channel the
     engine's closures sit in. The selection is an identification,
     made with no measured G.
  L6 (the third factor). Without the second moment as a factor, the
     lattice has 108 probability response readings. The whole of each
     list is one, with or without that factor.
  L7 (depth four). The engine's common mode is this whole on the
     searchlight's circulant. On the depth four web the engine reads
     this whole as the count of the lepton circulant, G M_b^2 =
     3 Q_l - 1 = 1 - 7.50 (alpha/2 pi)^2 (root.py, THE COMMON MODE;
     common_mode_count.py). The carrier is a node of the web and not an
     entry of the three lists.

RESULT. The coupling definition leaves c on a lattice, bounded below
by one, with no continuum. On this lattice the gates pass 71 values
with no measured G. Charge conjugation, the one nontrivial symmetry of
the modular data, keeps all four values that pass every convention.
The selection comes from fusion and universality. The closures the
common mode reads form a fusion closed set. Step 3 (every closure
reads the one common mode) requires that set to contain every channel
the engine's closures sit in. Only the whole does, so c = 1. This is
an identification made with no measured G. The boundary reading is the
coupling definition applied to a set of closures every one of which
the common mode reads, and G tests it.

Scope: the three lists of P2 at the boundary. The statement that the
common mode's coefficient is G M_b^2 is the framework's account of
what gravity is. No exchange dynamics at the boundary is computed.

Usage: python3 tests/probes/boundary_lattice.py
"""

import itertools
import math
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import boundary_whole as bw  # noqa: E402
from interference import root as r  # noqa: E402
from interference import mixing, wiring  # noqa: E402

G_SIGMA = r.G_REL_SIGMA                      # relative uncertainty of the recommended G


def lists():
    """The three lists of P2, each entry with its weight."""
    d = wiring.su3_qdims(3)
    D2 = round(sum(v * v for v in d.values()))
    s = bw.shares()
    return dict(
        probability={p: Fraction(round(d[p] ** 2), D2) for p in sorted(d)},
        response={"local": s["response"]["local"], "remainder": s["response"]["remainder"]},
        moment={"q=1": s["moment"]["mode"], "q=2": s["moment"]["whole"] - s["moment"]["mode"]})


def subsets(entries):
    """Every nonempty subset of a list, as a frozenset of its labels."""
    names = sorted(entries)
    return [frozenset(c) for n in range(1, len(names) + 1) for c in itertools.combinations(names, n)]


def weight(subset, entries):
    return sum(entries[x] for x in subset)


def run(report=print):
    L = lists()
    P, R, M = L["probability"], L["response"], L["moment"]
    assert sum(P.values()) == 1 and abs(sum(R.values()) - 1) < 1e-15 and abs(sum(M.values()) - 1) < 1e-12

    # L1: the lattice.
    p_sub, r_sub, m_sub = subsets(P), subsets(R), subsets(M)
    p_vals = sorted({weight(s, P) for s in p_sub})
    r_vals = sorted({weight(s, R) for s in r_sub})
    m_vals = sorted({round(weight(s, M), 12) for s in m_sub})
    cs = sorted({round(1 / (float(p) * rv * mv), 10) for p in p_vals for rv in r_vals for mv in m_vals})
    n_readings = len(p_vals) * len(r_vals) * len(m_vals)
    assert len(p_sub) == 1023 and p_vals == [Fraction(k, 36) for k in range(1, 37)]
    assert (len(r_vals), len(m_vals), n_readings) == (3, 2, 216)
    assert cs[0] == 1.0 and sum(c < 1 + 1e-9 for c in cs) == 1
    report(f"L1 {len(p_sub)} subsets of the primaries, {len(p_vals)} probabilities k/36, {n_readings} readings, "
           f"least c = {cs[0]:.6f}, attained once")

    # L2: no continuum.
    nearest = cs[1]
    shift = 1 / nearest - 1
    assert abs(nearest - 36 / 35) < 1e-9
    report(f"L2 the next reading is c = {nearest:.6f}: G {100 * shift:+.2f}%, {shift / G_SIGMA:+.0f} standard deviations")

    # L3: the gates with no measured G.
    wired = [c for c in cs if all(bw.gates(c))]
    six = [c for c in cs if bw.conventions_passed(c) == 6]
    n_wired = len(wired)
    assert (len(cs), n_wired) == (162, 71) and [round(c, 6) for c in six] == [1.0, round(36 / 35, 6), round(18 / 17, 6), round(12 / 11, 6)]
    report(f"L3 gates S and M as wired pass {n_wired} of {len(cs)} values. All six conventions of S pass c = "
           + ", ".join(f"{c:.4f}" for c in six))

    # L4: the symmetries of the modular data do not select.
    wzw0 = mixing._compute_wzw()
    S, idx = wzw0["S_mat"], wzw0["IDX"]
    conj = lambda p: (p[1], p[0])
    pairs = [(p, conj(p)) for p in sorted(P) if p < conj(p)]
    auts = []
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        perm = list(range(len(P)))
        for b, (a, c) in zip(bits, pairs):
            if b:
                perm[idx[a]], perm[idx[c]] = idx[c], idx[a]
        Pm = np.eye(len(P))[perm]
        if np.allclose(Pm @ S @ Pm.T, S):
            auts.append(bits)
    sym_values = sorted({round(1 / float(weight(x, P)), 9) for x in p_sub
                         if frozenset(conj(q) for q in x) == x})
    sym_six = [c for c in sym_values if bw.conventions_passed(c) == 6]
    assert auts == [(0,) * len(pairs), (1,) * len(pairs)]
    assert len(sym_values) == 22
    assert sym_six == [1.0, round(36 / 35, 9), round(18 / 17, 9), round(12 / 11, 9)]
    report(f"L4 symmetries of the modular data: {len(auts)} (identity, conjugation). Conjugation keeps "
           f"{len(sym_values)} values. All six conventions pass " + ", ".join(f"{c:.4f}" for c in sym_six))

    # L5: the fusion closed sets.
    wzw = mixing._compute_wzw()
    prim, N = wzw["PRIMARIES"], wzw["N_fus"]
    vac = wzw["IDX"][(0, 0)]
    closed = []
    for x in p_sub:
        idx = {wzw["IDX"][q] for q in x}
        if vac in idx and all(N[a, b, k] == 0 or k in idx for a in idx for b in idx for k in range(len(prim))):
            closed.append(x)
    fusion_c = sorted(1 / float(weight(x, P)) for x in closed)
    fund = (1, 0)
    holding = [x for x in closed if fund in x]
    assert [round(c, 9) for c in fusion_c] == [1.0, 3.0, 12.0, 36.0]
    assert len(holding) == 1 and holding[0] == frozenset(P)
    three = [x for x in closed if abs(1 / float(weight(x, P)) - 3) < 1e-9][0]
    assert all(bw.gates(3.0)) and not any(bw.gates(12.0)) and not all(bw.gates(36.0))
    assert all((p[0] - p[1]) % 3 == 0 for p in three)
    whole = float(weight(frozenset(P), P)) * sum(R.values()) * sum(M.values())
    assert abs(whole - 1) < 1e-12
    report(f"L5 fusion closed sets: c = " + ", ".join(f"{c:g}" for c in fusion_c)
           + ". The one holding the fundamental is the whole, G M_b^2 = 1. c = 3 passes the wired gates with no "
           "triality changer")

    # L5b: the engine's closures by channel, and gate U.
    h_of = lambda p: (p[0]**2 + p[1]**2 + p[0] * p[1] + 3 * p[0] + 3 * p[1]) / 3 / (r.k_SU3 + r.d11)
    assert abs(h_of((1, 0)) - float(r.h10)) < 1e-12 and abs(h_of((1, 1)) - float(r.h11)) < 1e-12
    channel = {"leptons (theta = h10)": float(r.h10), "(s,c,b) (h10/K^3)": float(r.h10),
               "(c,b,t) (h11/K^3)": float(r.h11)}
    reads = {}
    for x in closed:
        held = {name for name, h in channel.items() if any(abs(h_of(q) - h) < 1e-12 for q in x)}
        reads[round(1 / float(weight(x, P)), 9)] = sorted(held)
    universal = [c for c, held in reads.items() if len(held) == len(channel)]
    assert universal == [1.0] and reads[3.0] == ["(c,b,t) (h11/K^3)"]
    out_reads = reads
    report("L5b closures each fusion closed set reads: "
           + "; ".join(f"c = {c:g}: {', '.join(h) or 'none'}" for c, h in sorted(reads.items()))
           + ". Step 3 needs every closure read, and only the whole holds them all")

    # L6: without the second moment as a factor.
    assert len(p_vals) * len(r_vals) == 108
    report(f"L6 without the moment factor: {len(p_vals) * len(r_vals)} readings, the whole is one")

    # L7: the engine's depth four reading sits on this whole.
    frame, L = r.COMMON_MODE - whole, r.WEB.state["eq_loop"]
    searchlight = r.WEB["common_mode"].value_at(dict(r.WEB.state, eq_loop=0.0))
    assert abs(searchlight - whole) < 1e-14 and abs(frame / L + 7.502) < 0.001
    report(f"L7 the engine's common mode is this whole on the searchlight's circulant and the whole "
           f"{frame:+.4e} on the depth four web, the lepton circulant's count")

    report("RESULT: c is a lattice bounded below by one. The gates pass 71 values with no measured G, symmetry "
           "keeps all four that pass every convention, and fusion closure with every closure read keeps the whole, c = 1.")
    return dict(n_subsets=len(p_sub), n_readings=n_readings, n_values=len(cs), c_min=cs[0], nearest=nearest,
                nearest_sigma=shift / G_SIGMA, n_wired=n_wired, six=six, sym_values=sym_values,
                sym_six=sym_six, fusion_c=fusion_c, reads=out_reads, whole=whole, frame=frame)


if __name__ == "__main__":
    run()
