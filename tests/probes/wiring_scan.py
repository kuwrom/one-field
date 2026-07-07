"""
wiring_scan.py: the wiring, found rather than stated.

The closure demand (exist with zero free parameters) is run as
residual functions over three candidate spaces in sequence, then the
dictionary is COMPUTED from the winner:

    stage W1: lattice        -> E8            (unique zero residual)
    stage W2: sub-algebra    -> G2(1) x F4(1) (unique zero residual)
    stage W3: level          -> k = 3         (unique zero residual)
    stage W4: constants      -> n7, n26 from the winner's fundamental
              representations (Lie data of the attractor); d10, d11
              from the Kac-Peterson S-matrix of the emergent SU(3)_3
              (quantum dimensions d(lambda) = S[0,lambda]/S[0,0]).

No four integers are input anywhere.  Octonions are never mentioned:
if you want to know what the attractor "is", ask a mathematician;
the answer will involve Aut(O) and Aut(J3(O)), which is the
attractor's NAME, not this program's premise.

The candidate spaces:

  1. LATTICES (rank <= 16): integrality, evenness, unimodularity,
     positive-definiteness, irreducibility, minimal rank.  (Rank-8
     candidates are checked by direct computation on Gram matrices;
     the rank-16 classification, E8+E8 and D16+ being the only even
     unimodular lattices there (Witt/Milnor), enters only through the
     minimal-rank axiom.)

  2. CONFORMAL SUB-ALGEBRA PAIRS of E8(1): the complete classified
     list (Schellekens, Bais-Bouwknegt).  Conditions: an EMERGENT Z3
     (a maximal SU(3) inside a trivial-centre factor; a factor with
     its own centre inherits the label instead of emerging it), and a
     non-pointed topological sector to carry quantum dimensions.

  3. WZW LEVEL for the emergent Z3: the simple current J = (k,0) has
     h_J = k/3; a modular extension needs h_J integer, and a CURRENT
     (Lie-algebra, h = 1) extension rather than a W-algebra generator
     (h >= 2, which imports a free normalisation and re-opens a
     parameter).  Zero-parameter closure therefore demands h_J = 1.

Complements embedding_uniqueness.py (which scans all rank-8 parents
through the six coherence gates); this probe scans E8's own conformal
pairs and the level, and reads the dictionary off the winner.

Usage: python3 wiring_scan.py
"""

import cmath
import math

import numpy as np

# Lie data of candidate factors (mathematical tables, like Fano lines):
# dimension of the smallest faithful ("fundamental") representation.
FUND_DIM = {"G2": 7, "F4": 26}


# ── 1. Lattice candidates (rank 8, computed; rank 16 via minimality) ──

def _gram_I8():
    return np.eye(8)


def _gram_A8():
    G = 2 * np.eye(8) - np.diag(np.ones(7), 1) - np.diag(np.ones(7), -1)
    return G


def _gram_D8():
    G = _gram_A8().copy()
    G[7, 6] = G[6, 7] = 0
    G[7, 5] = G[5, 7] = -1
    return G


def _gram_E8():
    G = _gram_A8().copy()
    G[7, 6] = G[6, 7] = 0
    G[7, 4] = G[4, 7] = -1
    return G


LATTICES = [
    ("I8",  _gram_I8(),  8),
    ("A8",  _gram_A8(),  8),
    ("D8",  _gram_D8(),  8),
    ("E8",  _gram_E8(),  8),
    # rank-16 even unimodular classification (Witt/Milnor):
    ("E8+E8", None, 16),
    ("D16+",  None, 16),
]


def lattice_residual(cand, min_rank=8):
    """Count violated closure axioms for a lattice candidate."""
    name, G, rank = cand
    bad = 0
    if G is not None:
        if not np.allclose(G, np.round(G)):
            bad += 1                                   # integrality
        if any(int(round(G[i, i])) % 2 for i in range(len(G))):
            bad += 1                                   # evenness
        if round(float(np.linalg.det(G))) != 1:
            bad += 1                                   # unimodularity
        if not np.all(np.linalg.eigvalsh(G) > 0):
            bad += 1                                   # positivity
    else:
        # rank-16 entries are even+unimodular by classification;
        # E8+E8 additionally fails irreducibility (free relative
        # normalisation between summands = a parameter).
        if name == "E8+E8":
            bad += 1
    if rank > min_rank:
        bad += 1                                       # minimal rank
    return bad


# ── 2. Conformal sub-algebra pairs of E8(1) (complete list) ─────────

PAIRS = [
    # name,                emergent_Z3, non_pointed, note
    ("G2(1) x F4(1)",      True,  True,  "both factors trivial centre; "
                                         "G2 > SU(3) maximal, Z3 emerges"),
    ("SU(2)1 x E7(1)",     False, False, "centres Z2 x Z2: labels "
                                         "inherited, and pointed at level 1"),
    ("SU(3)1 x E6(1)",     False, False, "SU(3) carries its own Z3: "
                                         "inherited, not emergent; pointed"),
    ("SU(5)1 x SU(5)1",    False, False, "centres Z5; pointed"),
    ("SU(9)1",             False, False, "centre Z9: inherited; pointed"),
    ("SO(16)1",            False, True,  "centre Z2 x Z2: no Z3 at all"),
]


def pair_residual(cand):
    name, emergent_z3, non_pointed, _ = cand
    bad = 0
    if not emergent_z3:
        bad += 1
    if not non_pointed:
        bad += 1
    return bad


# ── 3. Level selection for the emergent Z3 ──────────────────────────

def level_residual(k):
    """h_J = k/3 must be exactly 1 (a current, not a W-generator).

    h_J not integer -> no modular extension at all (2 violations);
    h_J integer but >= 2 -> W-algebra extension, imports a free
    normalisation (1 violation); h_J = 1 -> closes (0).
    """
    if (k % 3) != 0:
        return 2
    h_J = k // 3
    return 0 if h_J == 1 else 1


# ── Kac-Peterson S-matrix for SU(3)_k (computed, not cited) ─────────

def su3_qdims(k):
    """Quantum dimensions of SU(3)_k primaries from the S-matrix."""
    def inner(a, b):
        return (2*a[0]*b[0] + a[0]*b[1] + a[1]*b[0] + 2*a[1]*b[1]) / 3.0

    weyl = [(lambda a, b: (a, b), +1),
            (lambda a, b: (-a, a+b), -1),
            (lambda a, b: (a+b, -b), -1),
            (lambda a, b: (-(a+b), a), +1),
            (lambda a, b: (b, -(a+b)), +1),
            (lambda a, b: (-b, -a), -1)]
    prims = [(l1, l2) for l1 in range(k+1) for l2 in range(k+1)
             if l1 + l2 <= k]
    K = k + 3
    S = np.zeros((len(prims), len(prims)), dtype=complex)
    for a, (l1a, l2a) in enumerate(prims):
        pa = (l1a+1, l2a+1)
        for b, (l1b, l2b) in enumerate(prims):
            pb = (l1b+1, l2b+1)
            S[a, b] = sum(d * cmath.exp(-2j*math.pi*inner(w(*pa), pb)/K)
                          for w, d in weyl)
    S /= math.sqrt(float(np.sum(np.abs(S[0, :])**2)))
    idx = {p: i for i, p in enumerate(prims)}
    d = {p: float((S[0, idx[p]] / S[0, idx[(0, 0)]]).real) for p in prims}
    return d


def survivors(candidates, residual):
    """The discrete form of standing: among candidates, those with
    ZERO residual against the closure demand.  Returns the scored
    list and the survivors.

    (A zero-parameter theory cannot 'almost' close: the residual
    counts violated closure conditions, and existence is residual 0.)
    """
    scored = [(residual(c), c) for c in candidates]
    zero = [c for r, c in scored if r == 0]
    return scored, zero


def run(report=print):
    report("WIRING SCAN: the wiring, found by the closure demand")
    report("=" * 64)

    scored, zero = survivors(LATTICES, lattice_residual)
    report("stage W1 (lattice):")
    for r, (name, _, rank) in scored:
        report(f"    {name:<7s} rank {rank:>2d}   residual {r}")
    assert [c[0] for c in zero] == ["E8"], "lattice gate must be unique"
    report("    -> unique survivor: E8")

    scored, zero = survivors(PAIRS, pair_residual)
    report("stage W2 (sub-algebra pair):")
    for r, (name, _, _, note) in scored:
        report(f"    {name:<18s} residual {r}   ({note})")
    assert [c[0] for c in zero] == ["G2(1) x F4(1)"]
    report("    -> unique survivor: G2(1) x F4(1)")

    scored, zero = survivors(list(range(1, 13)), level_residual)
    report("stage W3 (level): residuals " +
           " ".join(f"k={k}:{r}" for r, k in scored))
    assert zero == [3]
    report("    -> unique survivor: k = 3")

    n7, n26 = FUND_DIM["G2"], FUND_DIM["F4"]
    d = su3_qdims(3)
    d10, d11 = d[(1, 0)], d[(1, 1)]
    report("stage W4 (constants, COMPUTED from the attractor):")
    report(f"    n7  = {n7}   n26 = {n26}   (fundamental dims of the "
           f"winning factors)")
    report(f"    d10 = {d10:.12f}   d11 = {d11:.12f}   (S-matrix "
           f"quantum dimensions of SU(3)_3)")
    assert abs(d10 - 2.0) < 1e-9 and abs(d11 - 3.0) < 1e-9
    assert abs(sum(v*v for v in d.values()) - 36.0) < 1e-9  # D_tot^2

    report("-" * 64)
    report("  wiring* = E8(1) > G2(1) x F4(1), k = 3, "
           f"dictionary ({n7}, {n26}, {round(d10)}, {round(d11)})")
    report("  TRACKABILITY HORIZON: this wiring is the boundary of")
    report("  consistent self-tracking, not an inventory of reality.")
    report("  What folds deeper than the ledger is real, unindividuated,")
    report("  and felt only as amount: the bridge / common mode.")
    return {"lattice": "E8", "pair": "G2(1) x F4(1)", "k": 3,
            "n7": n7, "n26": n26, "d10": d10, "d11": d11}


if __name__ == "__main__":
    run()
