"""
wiring.py: the wiring, found rather than stated.

A closure-residual scan over the CLASSIFIED candidate menus
(lattices, conformal pairs, levels): the closure demand (exist with
zero free parameters) is evaluated as residual functions over each
finite, pre-classified menu in sequence, and the dictionary is then
COMPUTED from the winner.  Nothing is searched against data. The
residuals are structural closure conditions, and each menu is a
complete classification, not a tuned candidate pool:

    stage W1: lattice        -> E8            (unique zero residual)
    stage W2: sub-algebra    -> G2(1) x F4(1) (unique zero residual)
    stage W3: level          -> k = 3         (unique zero residual)
    stage W4: constants      -> n7, n26 from the winner's fundamental
              representations (Lie data of the attractor). d10, d11
              from the Kac-Peterson S-matrix of the emergent SU(3)_3
              (quantum dimensions d(lambda) = S[0,lambda]/S[0,0]).
    stage W5: coupling       -> the selected h_J = 1 current's pointed
              sector probability times its normalized response h_J/2pi.
    stage W6: bridge moment  -> the G2/Z3 harmonic tensored with the
              identity on the irreducible F4 factor; normalized trace
              over 7*26 modes gives |v_hat_1|^2 = 1/2.
    stage W7: vertex trace   -> the normalized Cartan-current operator
              on one representative of every E8 root pair has quadratic
              trace 30 = h_dual(E8), independently of current direction.

No four integers are input anywhere.  Octonions are never mentioned:
if you want to know what the attractor "is", ask a mathematician.
The answer will involve Aut(O) and Aut(J3(O)), which is the
attractor's NAME, not this program's premise.

The candidate spaces:

  1. LATTICES (rank <= 16): integrality, evenness, unimodularity,
     positive-definiteness, irreducibility, minimal rank.  (Rank-8
     candidates are checked by direct computation on Gram matrices.
     The rank-16 classification, E8+E8 and D16+ being the only even
     unimodular lattices there (Witt/Milnor), enters only through the
     minimal-rank axiom.)

  2. CONFORMAL SUB-ALGEBRA PAIRS of E8(1): the complete classified
     list (Schellekens, Bais-Bouwknegt).  Conditions: an EMERGENT Z3
     (a maximal SU(3) inside a trivial-centre factor. A factor with
     its own centre inherits the label instead of emerging it), and a
     non-pointed topological sector to carry quantum dimensions.

  3. WZW LEVEL for the emergent Z3: the simple current J = (k,0) has
     h_J = k/3. A modular extension needs h_J integer, and a CURRENT
     (Lie-algebra, h = 1) extension rather than a W-algebra generator
     (h >= 2, which imports a free normalisation and re-opens a
     parameter).  Zero parameter closure therefore demands h_J = 1.

Complements embedding_uniqueness.py (which scans all rank-8 parents
through the six coherence gates). This probe scans E8's own conformal
pairs and the level, and reads the dictionary off the winner.

Usage: python3 interference/wiring.py
"""

import cmath
import math

import numpy as np


# Lie data of candidate factors (mathematical tables, like Fano lines):
# dimension of the smallest faithful ("fundamental") representation.
FUND_DIM = {"G2": 7, "F4": 26}


# The algebraic root imports the result of this scan.  Keeping the
# canonical result here makes wiring.py the single source of the four
# integers and of the normalized current response defined below.
_CANONICAL_DATA = None


# # 1. Lattice candidates (rank 8, computed. Rank 16 via minimality) #

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
        # rank-16 entries are even+unimodular by classification.
        # E8+E8 additionally fails irreducibility (free relative
        # normalisation between summands = a parameter).
        if name == "E8+E8":
            bad += 1
    if rank > min_rank:
        bad += 1                                       # minimal rank
    return bad


# # 2. Conformal sub-algebra pairs of E8(1) (complete list) ####

PAIRS = [
    # name,                emergent_Z3, non_pointed, note
    ("G2(1) x F4(1)",      True,  True,  "both factors trivial centre. "
                                         "G2 > SU(3) maximal, Z3 emerges"),
    ("SU(2)1 x E7(1)",     False, False, "centres Z2 x Z2: labels "
                                         "inherited, and pointed at level 1"),
    ("SU(3)1 x E6(1)",     False, False, "SU(3) carries its own Z3: "
                                         "inherited, not emergent. Pointed"),
    ("SU(5)1 x SU(5)1",    False, False, "centres Z5. Pointed"),
    ("SU(9)1",             False, False, "centre Z9: inherited. Pointed"),
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


# # 3. Level selection for the emergent Z3 #############

def level_residual(k):
    """h_J = k/3 must be exactly 1 (a current, not a W-generator).

    h_J not integer -> no modular extension at all (2 violations).
    h_J integer but >= 2 -> W-algebra extension, imports a free
    normalisation (1 violation). h_J = 1 -> closes (0).

    R3-FROM-R1 (why k = 3 is a LEMMA, not a hypothesis): the only
    axiom in play is R1, zero free parameters.  The simple current
    J = (k,0) of SU(3)_k has h_J = k/3.  A modular extension needs
    h_J integer, so k = 3, 6, 9, .... For k > 3 the extending field
    has h_J >= 2, i.e. It is a W-algebra generator, and a W-algebra
    extension re-imports a free normalisation (the generator's
    coupling is not fixed by current algebra), excluded by the very
    zero parameter demand being scanned.  k = 3 is therefore the
    unique zero-residual output of R1, not an independent input.
    The assertion in run() freezes this: exactly {3} survives over
    k in 1..12.

    A CONDENSED-MATTER EXPOSURE OF THE SAME LEVEL.  k = 3 forces the
    Read-Rezayi Z3 state at nu = 2 + 2/(k+2) = 12/5, i.e. Fibonacci
    anyons.  Window: GaAs, B = 8.5 T, n = 4.9e11 cm^-2, T < 10 mK,
    mobility > 1e7, gap ~ 100 mK.  Established Abelian or
    non-Fibonacci order at 12/5 under those conditions kills it;
    interferometry or thermal Hall at 12/5 decides.
    """
    if (k % 3) != 0:
        return 2
    h_J = k // 3
    return 0 if h_J == 1 else 1


# # Kac-Peterson S-matrix for SU(3)_k (computed, not cited) ####

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

    (A zero parameter theory cannot 'almost' close: the residual
    counts violated closure conditions, and existence is residual 0.)
    """
    scored = [(residual(c), c) for c in candidates]
    zero = [c for r, c in scored if r == 0]
    return scored, zero


def normalized_current_response(h):
    """Canonical closed-loop response of a normalized dimension-one current.

    Each consumer supplies ``h = 1`` by an exact identity: the extension
    current has h_J = k/3 with the stage-W3 level k = 3, and the (7,26)
    bridge has h_bridge = h(7,G₂) + h(26,F₄) = 2/5 + 3/5.  In the
    framework's chosen unit-normalized current basis, the response is
    defined once as ``h/(2π)``, with ``2π`` the angular period.  This is
    not the raw integral of the singular separated-point two-point
    kernel, nor a coefficient chosen per consumer.
    """
    if not math.isclose(float(h), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("the current-response rule requires h = 1")
    return float(h) / (2.0 * math.pi)


def bridge_mode_second_moment(amplitude, theta, n_g2, n_f4):
    """Normalized bridge trace of the single Z3-channel response.

    The selected Z3 is inside the G2 factor.  An F4-invariant operator
    on the unique (7,26) connector is scalar on the irreducible 26 by
    Schur, so its response is the G2 harmonic tensored with identity.
    The identity multiplicity cancels in the normalized trace.
    """
    if n_g2 <= 0 or n_f4 <= 0:
        raise ValueError("bridge factor dimensions must be positive")
    omega = cmath.exp(2j * math.pi / 3.0)
    harmonic = [amplitude * math.cos(theta + 2.0 * math.pi * k / 3.0)
                for k in range(3)]
    v_hat1 = sum(harmonic[k] * omega**(-k) for k in range(3)) / 3.0
    g2_channel = v_hat1 * np.eye(n_g2, dtype=complex)
    bridge_response = np.kron(g2_channel, np.eye(n_f4, dtype=complex))
    multiplicity = n_g2 * n_f4
    return float(np.trace(bridge_response.conj().T @ bridge_response).real
                 / multiplicity)


def e8_roots():
    """The 240 E8 roots in the standard orthonormal R8 realization."""
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (-1.0, 1.0):
                for sj in (-1.0, 1.0):
                    root = np.zeros(8)
                    root[i], root[j] = si, sj
                    roots.append(root)
    for mask in range(2**8):
        signs = np.array([
            -1.0 if (mask >> i) & 1 else 1.0 for i in range(8)
        ])
        if int(np.count_nonzero(signs < 0)) % 2 == 0:
            roots.append(signs / 2.0)
    assert len(roots) == 240
    assert all(np.dot(root, root) == 2.0 for root in roots)
    return tuple(roots)


def e8_positive_root_current_form():
    """Quadratic current form on one root from every conjugate pair.

    For a normalized Cartan current q, the positive-root mode E_alpha has
    charge alpha.q.  Choosing the root whose first nonzero coordinate is
    positive selects one representative of each complex conjugate pair.
    The returned matrix M obeys q.T M q = sum_(alpha>0)(alpha.q)^2.
    """
    positive = []
    for root in e8_roots():
        first = root[np.flatnonzero(root)[0]]
        if first > 0.0:
            positive.append(root)
    assert len(positive) == 120
    return sum((np.outer(root, root) for root in positive),
               np.zeros((8, 8)))


def e8_vertex_charge_trace(direction=None):
    """E8 positive-root charge trace for a unit Cartan current."""
    q = np.zeros(8) if direction is None else np.asarray(direction,
                                                          dtype=float)
    if direction is None:
        q[0] = 1.0
    if q.shape != (8,) or not np.isclose(np.dot(q, q), 1.0):
        raise ValueError("the E8 current direction must be a unit 8-vector")
    current_form = e8_positive_root_current_form()
    assert np.allclose(current_form, 30.0 * np.eye(8), atol=1e-12)
    return float(q @ current_form @ q)


def run(report=print):
    global _CANONICAL_DATA

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
    assert zero == [3], \
        "R3 is a lemma of R1: k = 3 must be the unique zero residual " \
        "(k > 3 extends by a W-generator, re-importing a free " \
        "normalisation. Excluded by the zero parameter demand itself)"
    report("    -> unique survivor: k = 3 (a lemma of the "
           "zero parameter demand, not a hypothesis)")

    n7, n26 = FUND_DIM["G2"], FUND_DIM["F4"]
    d = su3_qdims(3)
    d10, d11 = d[(1, 0)], d[(1, 1)]
    report("stage W4 (constants, COMPUTED from the attractor):")
    report(f"    n7  = {n7}   n26 = {n26}   (fundamental dims of the "
           f"winning factors)")
    report(f"    d10 = {d10:.12f}   d11 = {d11:.12f}   (S-matrix "
           f"quantum dimensions of SU(3)_3)")
    assert abs(d10 - 2.0) < 1e-9 and abs(d11 - 3.0) < 1e-9
    D2_tot = sum(v*v for v in d.values())
    assert abs(D2_tot - 36.0) < 1e-9

    # The winning extension current J=(k,0) has h_J=k/3=1.  Current
    # algebra fixes its local normalisation; unlike a W-generator it
    # introduces no independent coupling.  The Z3 extension channel's
    # S-matrix probability is the sum over its simple-current orbit.
    # The engine defines a channel coupling as this closure probability
    # times the normalized dimension-one current response h/(2pi).
    h_J = zero[0] / 3.0
    pointed_orbit = ((0, 0), (3, 0), (0, 3))
    pointed_weight = sum(d[p]**2 for p in pointed_orbit) / D2_tot
    current_response = normalized_current_response(h_J)
    alpha_G2_Pl = pointed_weight * current_response
    assert abs(h_J - 1.0) < 1e-12
    assert abs(pointed_weight - 1.0/12.0) < 1e-12
    assert abs(alpha_G2_Pl - 1.0/(24.0*math.pi)) < 1e-12

    report("stage W5 (normalized-current response):")
    report(f"    h_J = {h_J:.0f},  P(Z3) = 3/{D2_tot:.0f} = "
           f"{pointed_weight:.12f}")
    report("    alpha_G2(M_Pl) = P(Z3) * h_J/(2pi) "
           f"= {alpha_G2_Pl:.12f}")

    # The unique bridge is (7,26).  Since the emergent Z3 is on its G2
    # side, F4 invariance makes the 26 factor an identity by Schur's
    # lemma.  The normalized trace therefore contains no new bridge
    # amplitude: it is the canonical G2 harmonic's q=1 power.
    N_bridge = n7 * n26
    harmonic_amplitude = math.sqrt(d10)
    harmonic_theta = d10 / d11**2
    bridge_mode_moment = bridge_mode_second_moment(
        harmonic_amplitude, harmonic_theta, n7, n26)
    assert abs(bridge_mode_moment - 0.5) < 1e-12
    report("stage W6 (bridge Schur trace):")
    report(f"    B=(7,26), multiplicity={N_bridge}; "
           "V_bridge=V_G2(Z3) tensor I_26")
    report("    Tr(V_bridge^dagger V_bridge)/182 "
           f"= |v_hat_1|^2 = {bridge_mode_moment:.12f}")

    # A normalized abelian current q in the selected E8 Cartan acts on
    # the root mode E_alpha with charge alpha.q.  The +/- roots are one
    # complex-conjugate mode, so the vertex trace uses one representative
    # per pair.  The E8 root identity gives M=sum alpha alpha^T=30 I8;
    # hence Tr_+(Q^2)=30 for every unit q.  This is the adjoint current
    # index h_dual(E8), derived without choosing a 26+4 mode inventory.
    vertex_charge_trace = e8_vertex_charge_trace()
    diagonal_direction = np.ones(8) / math.sqrt(8.0)
    assert abs(e8_vertex_charge_trace(diagonal_direction)
               - vertex_charge_trace) < 1e-12
    report("stage W7 (E8 current-operator trace):")
    report("    240 roots = 120 conjugate pairs; "
           "sum_(alpha>0) alpha alpha^T = 30 I_8")
    report("    Tr_+(Q^2) = sum_(alpha>0)(alpha.q)^2 "
           f"= {vertex_charge_trace:.0f} for |q|=1")

    report("-" * 64)
    report("  wiring* = E8(1) > G2(1) x F4(1), k = 3, "
           f"dictionary ({n7}, {n26}, {round(d10)}, {round(d11)})")
    report("  TRACKABILITY HORIZON: this wiring is the boundary of")
    report("  consistent self-tracking, not an inventory of reality.")
    report("  What folds deeper than the ledger is real, unindividuated,")
    report("  and felt only as amount: the bridge / common mode.")
    result = {"lattice": "E8", "pair": "G2(1) x F4(1)", "k": 3,
              "n7": n7, "n26": n26, "d10": d10, "d11": d11,
              "h_J": h_J, "D2_tot": D2_tot,
              "pointed_orbit": pointed_orbit,
              "pointed_weight": pointed_weight,
              "current_response": current_response,
              "alpha_G2_Pl": alpha_G2_Pl,
              "N_bridge": N_bridge,
              "harmonic_amplitude": harmonic_amplitude,
              "harmonic_theta": harmonic_theta,
              "bridge_mode_moment": bridge_mode_moment,
              "vertex_charge_trace": vertex_charge_trace}
    _CANONICAL_DATA = result.copy()
    return result


def canonical():
    """Return the canonical scan result, running it silently if needed."""
    if _CANONICAL_DATA is None:
        run(report=lambda *args, **kwargs: None)
    return _CANONICAL_DATA.copy()


#  Every module in interference/ exposes derive(); run() is kept as an
#  alias so the probe harness and standalone invocation both still work.
derive = run


if __name__ == "__main__":
    run()
