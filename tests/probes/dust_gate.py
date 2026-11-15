"""dust_gate.py: can the carrier itself supply the amount's dust?

THE QUESTION. dark_sector.py compares Omega_DM/Omega_b = 2 pi - 1 with
omega_c/omega_b, the amplitude of a component pressureless at
recombination (THE DUST GATE).  Dust is a conserved density carried
along a congruence with no pressure.  Its classical form in a field
theory is an identically conserved current built from gradients, the
comoving label current eps d phi^1 d phi^2 d phi^3 of the fluid
effective theory (Endlich et al. arXiv:1011.6396) and of mimetic
matter.  The stress of the traceless form itself is pure trace
(lambda_lane.py L6), so dust would have to come from such a current.
This probe asks whether the carrier, three complex components with
no names, can carry one.

PREMISES.  Field space C^3 = R^6 with its complex structure. The
amount carries no names, so a current that reads it is invariant under
relabelling the components (U(3) on the gradient bilinear's
components, all permutations included). The vacuum is the empty field,
translation invariant and isotropic, with no ruler profile (gate C,
lambda_scope.py S3).

STATEMENTS (the numbers below check each one).

  D1 (no name free constant 3-form).  The 3-forms on R^6 invariant
     under o(6) and under u(3) are zero. Under su(3) they are the two
     real parts of the holomorphic volume form, which carry charge 3
     under the overall phase and change sign under an odd relabelling
     of the components.

  D2 (every identically conserved current is a superpotential).  A
     current conserved without the field equations is the pullback of
     a closed 3-form on field space, J = *phi^*omega.  Field space is
     contractible, so omega = d beta and J^mu = d_nu S^{mu nu} with S
     antisymmetric and local (a trivial conservation law of the second
     kind: Olver, Applications of Lie Groups to Differential Equations,
     pp. 264-265).  Name free examples exist, e.g. from
     omega = d|phi|^2 ^ omega_K: J^mu = eps^{mu nu rho sigma}
     d_nu |phi|^2 d_rho J_sigma with J the common phase current, and
     S^{mu nu} = eps^{mu nu rho sigma} |phi|^2 d_rho J_sigma.  The
     superpotential of this example is checked by finite differences.
     That is a consistency check only, since d(rho dJ) is exact for any
     rho and J. The general statement is the cited theorem (Olver).

  D3 (no mean density).  J^0 = d_i S^{0i} is a spatial divergence, so
     its mean vanishes in every translation invariant state, at every
     order (lambda_scope.py S3). Its total charge is a boundary flux.
     A mean density needs a nonexact form. Such a form exists only on
     a noncontractible vacuum manifold, such as three locked phases
     (T^3, omega = d theta_1 ^ d theta_2 ^ d theta_3).  That current
     has density det(k) for windings theta_a = k_a . x. It therefore
     needs two things. The first is a winding profile, a ruler, which
     the empty field does not carry and which gives the graviton a mass
     (Carone arXiv:1911.13166). The second is an oriented triple of
     named components, a name (tests/probes/unnamed_metric.py).
     Computed: the density formula det(k). The ruler and name
     requirements are stated.

  D4 (the carrier's own charge is not the amount's dust).  The one
     current conserved on shell is the common phase charge.  On the gas
     sign its uniform state has p = epsilon, w = 1.  On the committed
     self-bound sign a uniform state at the zero pressure point has
     w = 0. But its sound speed is the phonon speed, c_s^2 = -a R0
     (0.757 c in the relativistic completion), five orders above the
     CMB bound. A space filling liquid is also not the empty field
     vacuum (tests/probes/self_bound_vacuum.py).  A dilute population of
     lumps is pressureless, but a lump is a closure with its clock. It
     belongs to the nameable side, not to the amount.  Stated.

RESULT.  The carrier supplies no dust for the amount. No name free
constant form exists. Every identically conserved current has zero
mean density in the empty field state. The only nonexact current needs
a ruler and a name. The conserved charge is either stiff (gas sign) or
a liquid whose sound speed is the phonon speed (committed sign).  The
dust gate cannot be passed at the classical level. Two routes remain.
One is the amount's induced stress, which needs its state. The other
is to read the amount as a conserved count that does no work, which
gives w = 0 and c_s = 0 as a declared premise.

Usage: python3 tests/probes/dust_gate.py     (a few seconds)
"""

import itertools

import numpy as np

EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


# ── D1: invariant 3-forms on R^6 = C^3 ────────────────────────────────

def _real_rep(A):
    """3x3 complex matrix -> 6x6 real matrix in the ordering (x1, x2, x3, y1, y2, y3)."""
    return np.block([[A.real, -A.imag], [A.imag, A.real]])


def _form_basis():
    triples = list(itertools.combinations(range(6), 3))
    basis = []
    for t in triples:
        w = np.zeros((6, 6, 6))
        for p in itertools.permutations(range(3)):
            w[tuple(t[i] for i in p)] = np.linalg.det(np.eye(3)[list(p)])
        basis.append(w)
    return basis


def _act(X, w):
    """Lie algebra action of a 6x6 generator on a 3-form."""
    return -(np.einsum("da,dbc->abc", X, w) + np.einsum("db,adc->abc", X, w)
             + np.einsum("dc,abd->abc", X, w))


def _invariants(generators):
    basis = _form_basis()
    M = np.array([np.concatenate([_act(X, w).ravel() for X in generators]) for w in basis]).T
    _, s, vt = np.linalg.svd(M)
    null = vt[np.sum(s > 1e-9):]
    return [sum(c * w for c, w in zip(v, basis)) for v in null]


def d1():
    o6 = []
    for i, j in itertools.combinations(range(6), 2):
        X = np.zeros((6, 6)); X[i, j], X[j, i] = 1.0, -1.0
        o6.append(X)
    su3 = []
    for i, j in itertools.combinations(range(3), 2):
        for z in (1.0, 1j):
            A = np.zeros((3, 3), complex); A[i, j], A[j, i] = z, -np.conj(z)
            su3.append(_real_rep(A))
    for k in range(2):
        A = np.zeros((3, 3), complex); A[k, k], A[k + 1, k + 1] = 1j, -1j
        su3.append(_real_rep(A))
    phase = _real_rep(1j * np.eye(3, dtype=complex))
    inv_o6, inv_u3, inv_su3 = _invariants(o6), _invariants(su3 + [phase]), _invariants(su3)
    # the su(3) invariants rotate into each other under the overall phase, with charge 3
    a, b = inv_su3
    Pa = _act(phase, a)
    coeff = [np.sum(Pa * a) / np.sum(a * a), np.sum(Pa * b) / np.sum(b * b)]
    charge = float(np.hypot(*coeff))
    swap = np.eye(6)[[1, 0, 2, 4, 3, 5]]                           # relabel components 1 <-> 2
    odd = [float(np.sum(np.einsum("ad,be,cf,def->abc", swap, swap, swap, w) * w) / np.sum(w * w))
           for w in inv_su3]
    return dict(o6=len(inv_o6), u3=len(inv_u3), su3=len(inv_su3), phase_charge=charge, odd_relabel=odd)


# ── D2: a name free identically conserved current and its superpotential ──

def _field(rng, n_waves=3):
    k = rng.normal(size=(3, n_waves, 4))
    a = rng.normal(size=(3, n_waves)) + 1j * rng.normal(size=(3, n_waves))
    return lambda x: np.array([np.sum(a[c] * np.exp(1j * (k[c] @ x))) for c in range(3)])


def _grad(f, x, h=1e-3):
    out = []
    for mu in range(4):
        e = np.zeros(4); e[mu] = h
        out.append((-f(x + 2 * e) + 8 * f(x + e) - 8 * f(x - e) + f(x - 2 * e)) / (12 * h))
    return np.array(out)


def d2(rng, n_points=5):
    phi = _field(rng)
    rho = lambda x: float(np.sum(np.abs(phi(x))**2))
    Jc = lambda x: np.imag(np.conj(phi(x)) @ _grad(phi, x).T)                  # common phase current
    dJ = lambda x: _grad(Jc, x)                                               # dJ[rho, sigma] = d_rho J_sigma
    cur = lambda x: np.einsum("mnrs,n,rs->m", EPS4, _grad(rho, x), dJ(x))      # eps d rho dJ
    S = lambda x: np.einsum("mnrs,rs->mn", EPS4, rho(x) * dJ(x))              # eps rho dJ
    worst_div, worst_sup = 0.0, 0.0
    for _ in range(n_points):
        x = rng.normal(size=4)
        c = cur(x)
        div = np.trace(_grad(cur, x))
        sup = np.einsum("nmn->m", _grad(S, x))                                 # d_nu S^{mu nu}
        scale = max(1.0, float(np.abs(c).max()))
        worst_div = max(worst_div, abs(div) / scale)
        worst_sup = max(worst_sup, float(np.abs(sup - c).max()) / scale)
    return dict(divergence=worst_div, superpotential=worst_sup)


# ── D3: the nonexact current needs a winding profile ─────────────────

def d3(rng):
    k = rng.normal(size=(3, 3))                                   # theta_a = k_a . x
    density = float(np.linalg.det(k))                             # eps grad th1 grad th2 grad th3
    return dict(winding_density=density, no_winding=float(np.linalg.det(np.zeros((3, 3)))))


def run(report=print, seed=11):
    rng = np.random.default_rng(seed)
    r1, r2, r3 = d1(), d2(rng), d3(rng)
    assert r1["o6"] == 0 and r1["u3"] == 0 and r1["su3"] == 2
    assert abs(r1["phase_charge"] - 3.0) < 1e-9
    assert all(abs(o + 1.0) < 1e-9 for o in r1["odd_relabel"])
    assert r2["divergence"] < 1e-5 and r2["superpotential"] < 1e-5
    assert abs(r3["winding_density"]) > 1e-3 and r3["no_winding"] == 0.0
    report(f"D1 invariant 3-forms on R^6: o(6) {r1['o6']}, u(3) {r1['u3']}, su(3) {r1['su3']} "
           f"(charge {r1['phase_charge']:.1f} under the overall phase, sign {r1['odd_relabel'][0]:+.0f} "
           f"under an odd relabelling): no name free constant form")
    report(f"D2 name free conserved current eps d|phi|^2 dJ: divergence {r2['divergence']:.1e}, "
           f"equals d_nu S^(mu nu) to {r2['superpotential']:.1e}: a superpotential, zero mean density "
           f"in a translation invariant state")
    report(f"D3 the nonexact current on T^3 has density det(k) = {r3['winding_density']:+.3f} for a "
           f"winding profile and {r3['no_winding']:.0f} without one: it needs a ruler and a named triple")
    report("D4 the common phase charge: w = 1 on the gas sign. On the committed sign w = 0 at the "
           "zero pressure point but c_s^2 = -a R0 (0.757 c relativistic) against c_s^2 < 3.2e-6, and it "
           "fills space. Lumps are pressureless but are closures, not the amount")
    report("RESULT: the carrier supplies no dust for the amount. The dust gate is not passed at the "
           "classical level. What remains is the amount's induced stress, or its reading as a "
           "conserved count that does no work (w = 0, c_s = 0) as a declared premise.")
    return dict(d1=r1, d2=r2, d3=r3)


if __name__ == "__main__":
    run()
