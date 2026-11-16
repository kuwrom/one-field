"""The stabilizer of the saddle: why every row reads one cone.

THE QUESTION. The summed kernel of ledger_kernel.py and
one_list_kernel.py expands every row about one flat metric, so one
metric for all rows is their premise P1, and an induced gravity needs
it: "all the linearized background fields must see the 'same'
Lorentzian metric", while "multiple Lorentzian metrics (and worse)
seems generic" (Barcelo, Liberati, Visser arXiv:gr-qc/0111059 Sec.
VII). With no background metric the
vacuum holds one index carrying datum, the constant saddle h_vac of
signature (1,3) (signature_of_the_saddle.py S2, S3). What can a row's
kinetic operator, a second composite metric or a regulator be built
from there, and do scalar, Weyl and vector rows share a cone?

PREMISES. P1 the vacuum is the empty field on the saddle
h_vac = diag(a0, -a1, -a1, -a1): it carries no other tensor, so its
symmetry is the stabilizer of h_vac in gl(4). P2 a row's quadratic
operator, a composite rank-2 tensor's vacuum value and a regulator's
weight are invariant under that symmetry.

STATEMENTS (the numbers below check each one).

  Z1 (the stabilizer). {X : X^T h + h X = 0} has dimension 6 for every
     nondegenerate saddle, so(1,3) on the passing one. Its invariant
     tensors: no vector, one rank-2 tensor (h_vac itself, symmetric; no
     antisymmetric one), none of rank 3, four of rank 4 (three products
     of h_vac and the Levi-Civita density).
  Z2 (scalar rows). The second order symbol K^{mu nu} p_mu p_nu has
     K = k h_vac^-1: every scalar row vanishes on h_vac^-1(p, p) = 0.
  Z3 (vector rows). The symbol K^{mu nu alpha beta} p_alpha p_beta is
     a p^2 h^-1 + b p p, the Levi-Civita term dropping against p p:
     its determinant is a^3 (a + b) (p^2)^4 / det h_vac, zero only on
     the same cone.
  Z4 (Weyl rows). The first order symbol M^mu p_mu between the two
     spinor representations of the stabilizer has one invariant, and
     its determinant is c p^2: the same cone. Z2 to Z4 are the
     condition of mono metricity, "the Fresnel equation must be some
     integer power of some single quadratic expression" (Barcelo,
     Liberati, Visser arXiv:gr-qc/0505065 Sec. 7.2).
  Z5 (a second composite metric). Any symmetric rank-2 vacuum value is
     a multiple of h_vac, so a metric built from other fields (a spinor
     bilinear's vierbein, Wetterich hep-th/0307145) has the cone of
     h_vac. The multiple is a field normalization.
  Z6 (the regulator). A weight invariant under the stabilizer is a
     function of h_vac^-1(p, p) alone and holds no frame vector (Z1),
     so the frame part D of signature_of_the_saddle.py S4 vanishes.
  Z7 (what a frame would do). Add one vector u to the vacuum: the
     stabilizer falls to so(3), dimension 3, with one invariant vector
     and two symmetric rank-2 tensors. A scalar row's symbol is then
     a h^-1 + b u u and a Weyl row's has two invariants: each row
     carries its own speed, sqrt(a1/a0) no longer shared. This is the
     case of a uniform clock filling space (carrier_sakharov.py: order
     one in c_T, 5/7 under a spatial cutoff) and of a preferred frame
     at the cutoff, where "different fields have different values of
     c, with fractional differences around 0.1% to 10%" (Collins,
     Perez, Sudarsky, Urrutia, Vucetich arXiv:gr-qc/0403053). Running
     does not repair it for the amount: "the running of the limiting
     velocities only happens due to the interactions that couple one
     type of particle to another" (Anber and Donoghue
     arXiv:1102.0789), and the amount has none. Its cone is the
     common one by Z2 alone.
  Z8 (weight is energy). A row couples to the vacuum only through h, so
     its response to a change of h, the tensor T = 2 dS/dh that the kernel's
     pole couples to, is its whole kinetic form. For a wave of a scalar
     row, at rest and moving, the period average of T^0_0 equals the
     canonical energy density, and T^0_x the canonical momentum
     density: what gravitates is what resists. With the frame term of
     Z7 the ratio is a h^00/(a h^00 + b), below one and different for
     each row: an Eotvos signal of order b. The Weyl row, a wave of
     its one invariant symbol, and the vector row, a wave of Maxwell's
     form on the saddle, give the same: ratio one and the saddle's
     speed with no frame term, and with it a e^0/(a e^0 + b) and
     1 - b/(2 c_E), c_E the coefficient of the electric term.

RESULT. On the empty field's saddle the rows cannot disagree about the
cone: representation theory leaves each row's symbol one invariant, a
multiple of h_vac^-1(p, p) or a root of it, and the tensor the metric
couples to is then the row's energy and momentum. One cone and equal
weight are the statement that the vacuum holds only h_vac.
Both fail exactly when the vacuum, or a regulator, holds a vector.

Scope: constant saddle, quadratic operators with at most two
derivatives, the stabilizer's connected component. The counts are
exact nullspace dimensions (singular values below 1e-10).

Usage: python3 tests/probes/saddle_stabilizer.py     (under a second)
"""

import itertools

import numpy as np

TOL = 1e-10
SADDLE = np.diag([1.7, -0.6, -0.6, -0.6])          # a0, a1 arbitrary: nothing below reads their values
U = np.array([1.0, 0.0, 0.0, 0.0])                 # a frame vector, for Z7
SIGMA = [np.eye(2, dtype=complex),
         np.array([[0, 1], [1, 0]], dtype=complex),
         np.array([[0, -1j], [1j, 0]], dtype=complex),
         np.array([[1, 0], [0, -1]], dtype=complex)]
EPS = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS[_p] = np.linalg.det(np.eye(4)[list(_p)])


def nullspace(A):
    """Orthonormal basis of the kernel of A, as columns."""
    _, s, vh = np.linalg.svd(A)
    rank = int((s > TOL * max(1.0, s[0])).sum()) if s.size else 0
    return vh[rank:].conj().T


def stabilizer(h, vectors=()):
    """Basis of {X in gl(4): X^T h + h X = 0, X v = 0 for each fixed vector v}, X acting on vectors."""
    rows = []
    for a in range(16):
        X = np.zeros(16)
        X[a] = 1.0
        X = X.reshape(4, 4)
        rows.append(np.concatenate([(X.T @ h + h @ X).ravel()] + [X @ v for v in vectors]))
    ns = nullspace(np.array(rows).T)
    return [ns[:, i].reshape(4, 4) for i in range(ns.shape[1])]


def act(X, rank, upper):
    """The generator X on tensors of the given rank, indices all upper (vectors) or all lower (covectors)."""
    G = X if upper else -X.T
    total = np.zeros((4**rank, 4**rank))
    for slot in range(rank):
        factors = [np.eye(4)] * rank
        factors[slot] = G
        m = factors[0]
        for f in factors[1:]:
            m = np.kron(m, f)
        total += m
    return total


def invariants(algebra, rank, upper=True):
    """Basis (as flat columns) of the tensors of that rank annihilated by every generator."""
    return nullspace(np.vstack([act(X, rank, upper) for X in algebra]))


def spinor(X, h):
    """Left handed spinor generator of X in so(h), through a frame e with e^T diag(1,-1,-1,-1) e = h."""
    e = np.diag(np.sqrt(np.abs(np.diag(h))))               # diagonal saddle: frame components
    w = e @ X @ np.linalg.inv(e)                          # the generator in the orthonormal frame
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    wl = eta @ w                                          # omega_{ab}, antisymmetric
    sbar = [SIGMA[0]] + [-s for s in SIGMA[1:]]
    S = np.zeros((2, 2), dtype=complex)
    for a in range(4):
        for b in range(4):
            S += 0.25 * wl[a, b] * 0.5 * (SIGMA[a] @ sbar[b] - SIGMA[b] @ sbar[a])
    return S


def weyl_invariants(algebra, h):
    """Symbols M^mu (2x2 each) with psi^dagger M^mu p_mu psi invariant: columns of a (16)-dim complex basis."""
    blocks = []
    for X in algebra:
        S = spinor(X, h)
        # psi -> (1 + S) psi and p -> (1 - X^T) p leave psi^dagger M^mu p_mu psi fixed when
        # S^dagger M^mu + M^mu S - X^mu_nu M^nu = 0, on M flattened as (mu, i, j)
        op = (np.kron(np.eye(4), np.kron(S.conj().T, np.eye(2)))        # S^dagger on the row index
              + np.kron(np.eye(4), np.kron(np.eye(2), S.T))             # S on the column index
              - np.kron(X, np.eye(4)))
        blocks.append(op)
    return nullspace(np.vstack(blocks))


def weight_over_energy(a, b, hinv, k, m=1.0, n=4096):
    """Period averages for the wave cos(w t - k x) of L = (1/2) K^{mu nu} d phi d phi - (1/2) m^2 phi^2.

    K = a h^-1 + b u u with u = dt.  Returns (T^0_0 / canonical energy,
    T^0_x / canonical momentum): T = 2 dS/dh varies h alone, the frame
    term b u u does not move with it.
    """
    Ktt, Kxx = a * hinv[0, 0] + b, a * hinv[1, 1]
    w = np.sqrt((m * m - Kxx * k * k) / Ktt)                    # K(p, p) = m^2 with p = (w, -k)
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    phi, dt, dx = np.cos(th), -w * np.sin(th), k * np.sin(th)       # th = w t - k x
    L = 0.5 * (Ktt * dt**2 + Kxx * dx**2) - 0.5 * m * m * phi**2
    energy = (Ktt * dt * dt - L).mean()                          # pi phi_dot - L
    momentum = (Ktt * dt * dx).mean()                            # pi d_x phi
    weight = (a * hinv[0, 0] * dt * dt - L).mean()               # T^0_0 of the h-built part
    flow = (a * hinv[0, 0] * dt * dx).mean()                     # T^0_x
    return weight / energy, (flow / momentum if k else 1.0)


def weyl_wave(a, b, hinv, kvec):
    """A Weyl wave of M^0 = (a e0 + b) 1, M^i = a e1 sigma^i: (speed, T^0_0 / canonical energy).

    e0, e1 are the frame's components on the diagonal saddle. The frame
    term b does not move with h, so T^0_0 = a e0 psi^dagger psi w.
    """
    e0, e1 = np.sqrt(hinv[0, 0]), np.sqrt(-hinv[1, 1])
    sk = sum(k * s for k, s in zip(kvec, SIGMA[1:]))
    w, v = np.linalg.eigh(a * e1 * sk / (a * e0 + b))            # (M^0 w - a e1 sigma.k) psi = 0
    psi, omega = v[:, -1], w[-1]
    norm = float(np.real(psi.conj() @ psi))
    return omega / np.linalg.norm(kvec), (a * e0 * norm * omega) / ((a * e0 + b) * norm * omega)


def vector_wave(a, b, hinv, k=0.7, n=4096):
    """A transverse Maxwell wave on the saddle with a frame term b E^2/2: (speed, T^0_0 / canonical energy)."""
    cE, cB = a * hinv[0, 0] * -hinv[1, 1] + b, a * hinv[1, 1]**2
    w = k * np.sqrt(cB / cE)
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    E, Bf = w * np.sin(th), k * np.sin(th)                      # A = cos(w t - k x) along y
    energy = (0.5 * (cE * E * E + cB * Bf * Bf)).mean()
    weight = (0.5 * ((cE - b) * E * E + cB * Bf * Bf)).mean()
    return w / k, weight / energy


def run(report=print):
    out = {}
    rng = np.random.default_rng(7)
    h, hinv = SADDLE, np.linalg.inv(SADDLE)
    p2 = lambda p: float(p @ hinv @ p)
    null = np.array([np.sqrt(-hinv[1, 1] / hinv[0, 0]) * np.linalg.norm(k := rng.normal(size=3)), *k])
    generic = rng.normal(size=4)
    assert abs(p2(null)) < 1e-12 < abs(p2(generic))

    # Z1: the stabilizer and its invariant tensors.
    g = stabilizer(h)
    dims = {sig: len(stabilizer(np.diag(sig))) for sig in itertools.product((1.0, -1.0), repeat=4)}
    count = {(k, up): invariants(g, k, up).shape[1] for k in (1, 2, 3, 4) for up in (True, False)}
    inv2 = invariants(g, 2, False)[:, 0].reshape(4, 4)
    assert len(g) == 6 and set(dims.values()) == {6}
    assert count == {(1, True): 0, (1, False): 0, (2, True): 1, (2, False): 1, (3, True): 0, (3, False): 0,
                     (4, True): 4, (4, False): 4}
    assert np.abs(inv2 - inv2.T).max() < TOL and np.abs(inv2 / inv2[0, 0] - h / h[0, 0]).max() < 1e-9
    out.update(stabilizer_dim=len(g), invariant_counts={f"rank {k}": count[(k, False)] for k in (1, 2, 3, 4)})
    report(f"Z1 stabilizer of the saddle: dimension {len(g)} (6 for all 16 sign patterns). Invariant tensors of "
           f"rank 1/2/3/4: {count[(1, False)]}/{count[(2, False)]}/{count[(3, False)]}/{count[(4, False)]}. The "
           "rank-2 one is the saddle, symmetric")

    # Z2: scalar rows.
    K = invariants(g, 2, True)[:, 0].reshape(4, 4)
    assert np.abs(K / K[0, 0] - hinv / hinv[0, 0]).max() < 1e-9
    report("Z2 scalar rows: the second order symbol is a multiple of the inverse saddle")

    # Z3: vector rows, a random invariant rank-4 symbol contracted with p p.
    basis4 = invariants(g, 4, True)
    def rel_det(T, p):
        sym = np.einsum("manb,a,b->mn", T, p, p)
        return abs(np.linalg.det(sym)) / np.abs(sym).max() ** 4

    symbols = [(basis4 @ rng.normal(size=basis4.shape[1])).reshape(4, 4, 4, 4) for _ in range(20)]
    dets = {"null": max(rel_det(T, null) for T in symbols),
            "generic": float(np.median([rel_det(T, generic) for T in symbols]))}
    eps_term = np.einsum("manb,a,b->mn", EPS, generic, generic)
    assert basis4.shape[1] == 4 and dets["null"] < 1e-12 and dets["generic"] > 1e-3
    assert np.abs(eps_term).max() < TOL
    out.update(vector_det=dets)
    report(f"Z3 vector rows: 4 invariant symbols, the Levi-Civita one drops against p p. Determinant on the cone "
           f"at most {dets['null']:.1e}, off it {dets['generic']:.1e} at the median of 20 symbols (relative)")

    # Z4: Weyl rows.
    W = weyl_invariants(g, h)
    assert W.shape[1] == 1
    M = W[:, 0].reshape(4, 2, 2)
    det_w = lambda p: np.linalg.det(np.einsum("mij,m->ij", M, p))
    ratio = det_w(generic) / p2(generic)
    other = rng.normal(size=4)
    assert abs(det_w(null)) < 1e-12 and abs(det_w(other) / p2(other) - ratio) < 1e-9 * abs(ratio) > 0
    out.update(weyl_invariants=W.shape[1])
    report(f"Z4 Weyl rows: {W.shape[1]} invariant first order symbol. Its determinant is a constant times "
           "h^-1(p, p)")

    # Z5, Z6: a second composite metric and the regulator read the same invariants.
    assert count[(2, False)] == 1 and count[(1, True)] == 0
    report("Z5 any symmetric rank-2 vacuum value is a multiple of the saddle: a second composite metric shares "
           "its cone.  Z6 no invariant vector: a symmetric regulator holds no frame part")

    # Z7: a frame vector in the vacuum.
    gf = stabilizer(h, vectors=(U,))
    cf = {k: invariants(gf, k, True).shape[1] for k in (1, 2)}
    basis2 = invariants(gf, 2, True)
    sym_dim = np.linalg.matrix_rank(np.array([(v.reshape(4, 4) + v.reshape(4, 4).T).ravel() for v in basis2.T]),
                                    tol=1e-8)
    Wf = weyl_invariants(gf, h)
    speeds = []
    for _ in range(4):                                     # scalar rows a h^-1 + b u u: each its own speed
        a, b = rng.uniform(0.5, 2.0, size=2)
        Kf = a * hinv + b * np.outer(U, U)
        speeds.append(float(np.sqrt(-Kf[1, 1] / Kf[0, 0])))
    shared = float(np.sqrt(-hinv[1, 1] / hinv[0, 0]))
    assert len(gf) == 3 and cf[1] == 1 and sym_dim == 2 and Wf.shape[1] == 2
    assert max(speeds) - min(speeds) > 1e-2 and all(abs(s - shared) > 1e-3 for s in speeds)
    out.update(frame_stabilizer_dim=len(gf), frame_vectors=cf[1], frame_symmetric_rank2=int(sym_dim),
               frame_weyl_invariants=Wf.shape[1], frame_speeds=speeds, shared_speed=shared)
    report(f"Z7 with a frame vector: stabilizer dimension {len(gf)}, {cf[1]} invariant vector, {sym_dim} symmetric "
           f"rank-2 tensors, {Wf.shape[1]} Weyl symbols. Four scalar rows then run at "
           + ", ".join(f"{s:.3f}" for s in speeds) + f" against the saddle's {shared:.3f}")

    # Z8: the tensor the metric couples to against the canonical energy and momentum.
    rows = [(float(a), float(b)) for a, b in rng.uniform(0.5, 2.0, size=(4, 2))]
    exact = [weight_over_energy(a, 0.0, hinv, k) for a, _ in rows for k in (0.0, 0.3, 0.6)]
    framed = [weight_over_energy(a, b, hinv, 0.0)[0] for a, b in rows]
    expect = [a * hinv[0, 0] / (a * hinv[0, 0] + b) for a, b in rows]
    assert max(abs(x - 1.0) for pair in exact for x in pair) < 1e-12
    assert max(abs(f - e) for f, e in zip(framed, expect)) < 1e-12
    eotvos = 2 * (max(framed) - min(framed)) / (max(framed) + min(framed))
    assert max(framed) < 1.0 and eotvos > 1e-2
    kvec = rng.normal(size=3)
    weyl = [weyl_wave(a, 0.0, hinv, kvec) for a, _ in rows]
    weyl_f = [weyl_wave(a, b, hinv, kvec) for a, b in rows]
    vec = [vector_wave(a, 0.0, hinv) for a, _ in rows]
    vec_f = [vector_wave(a, b, hinv) for a, b in rows]
    other = max(max(abs(s - shared), abs(x - 1.0)) for s, x in weyl + vec)
    assert other < 1e-12
    assert all(x < 0.99 and abs(s - shared) > 1e-2 for s, x in weyl_f + vec_f)
    assert max(x for _, x in weyl_f) - min(x for _, x in weyl_f) > 1e-2
    out.update(weight_over_energy=max(abs(x - 1.0) for pair in exact for x in pair), frame_weight=framed,
               frame_eotvos=eotvos, weyl_vector_error=other, frame_weight_weyl=[x for _, x in weyl_f],
               frame_weight_vector=[x for _, x in vec_f])
    report(f"Z8 weight over energy and flow over momentum: 1 to {out['weight_over_energy']:.0e} for four scalar "
           "rows at rest and moving. With the frame term "
           + ", ".join(f"{f:.3f}" for f in framed) + f", an Eotvos parameter {eotvos:.2f} between rows; Weyl and "
           f"vector waves: speed and weight at the saddle's values to {other:.0e}, with the frame term weights "
           + ", ".join(f"{x:.3f}" for _, x in weyl_f) + " and " + ", ".join(f"{x:.3f}" for _, x in vec_f))
    report("RESULT: on the empty field's saddle every row's symbol is one invariant, so the rows share the cone "
           "and their weight is their energy; a vector in the vacuum or in a regulator gives each row "
           "its own speed and its own weight.")
    return out


if __name__ == "__main__":
    run()
