"""The vacuum term, derived: the traceless form is the only lane that closes.

THE AXIOM. To be something is to carry non zero interference weight.
A constant c added to the potential, the vacuum term, changes no
relative phase between configurations of the field.  Interference is
only about relative phase, so the constant carries no interference
weight.  By the axiom it must then make no difference at all.  In
particular it may not weigh one history against another, so its
contribution to the phase must be the same for every history.  This is
the condition of history independence.

THE CLASSIFICATION. Let U enter a local action, algebraic in U and
first order in the fields' derivatives, through L(h, U).  The shift
U -> U + c adds c dL/dU + O(c^2).  History independence for every c
needs d^2L/dU^2 = 0 and dL/dU = -mu with c int mu the same for every
history.  Exactly three densities do that:
  (i)   mu = 0: the action carries no potential;
  (ii)  mu a closed four form built from the fields, the pullback of a
        closed form Omega on field space: then c int mu is a boundary
        term (the two measure and Pontryagin constructions);
  (iii) mu a fixed density: then c int mu = c Vol, the same for every
        history.
The full form sqrt|det h| / U has d^2L/dU^2 = 2 sqrt|det h| / U^3 != 0.
There the constant weighs each history by its own measure. That is how
Lambda gravitates, and it means the full form fails the axiom before
any data are compared.

STATEMENTS (the numbers below check each one).

  L1 (the forms against the axiom). Across random configurations,
     d^2L/dU^2 is nonzero for the full form. For the traceless form
     dL/dU = -1 identically, and for the two measure form dL/dU = -Phi.
  L2 (case (ii) closes, for every Omega). Take a random constant closed
     4-form Omega on the six real components. Wherever the pullback is
     nonzero, the gradients span all four directions, so h is positive
     definite (Euclidean).  Take a Gaussian state with a Lorentz
     invariant correlator, <d_mu phi^i d_nu phi^j> = S^ij eta_mu nu
     (the only invariant two index tensor is eta). Wick's theorem then
     contracts the antisymmetric Omega with the symmetric S, so
     <mu> = 0.  A nonzero value needs an antisymmetric spacetime part
     B_mu nu, which fixes a frame.  The Euclidean statement holds for
     every pullback, for any coefficients of Omega. The expectation
     value is checked for constant coefficients.  The rank argument does
     not cover the same case built from an independent gauge field, the
     Pontryagin density (arXiv:2301.01662 Eqs. 50-58). That case needs
     <F F~> != 0, a CP-odd condensate that a vacuum with theta bar = 0
     does not carry (interference/couplings.py).  So that closure is
     conditional on the strong CP reading, and the uniqueness below is
     conditional on it too.
  L3 (case (i) closes). With no potential the carrier has no mass term
     and no self-bound sign: no lump, no clock (lump_in_empty_field.py,
     self_bound_vacuum.py).  Stated.
  L4 (case (iii) is one form). With a fixed density the remaining
     freedom is f(det h). A kinetic term quadratic in the gradients
     needs f homogeneous of degree two in them. det h has degree eight,
     so f = c |det h|^(1/4). This is the traceless form
     2 |det h|^(1/4) - U, with the auxiliary unimodular metric
     eliminated (composite_metric.py C2).
  L2b (field dependent coefficients). For Omega(phi) the Wick pairings
     also join a field inside the coefficient to a derivative.  In a
     Lorentz invariant state that correlator is a vector with nothing to
     point along, so <phi^m d_mu phi^i> = 0. The derivative pairs are
     S eta again. So <mu> = 0 for every coefficient function.  Checked
     by the exact Isserlis sum over all 15 pairings, for coefficients
     quadratic in the fields. A state with a vector v_mu and an
     antisymmetric part gives a nonzero value.  Beyond Gaussian order:
     the action is built from h and U, both parity even, so in a parity
     symmetric vacuum the pseudoscalar <mu> vanishes at every order.
     Stated.
  L5 (the lane closes). About a constant saddle h_vac, the first order
     term of the traceless form is (1/2)|det h_vac|^(1/4) sum_mu
     (1/sigma_mu) |d_mu phi|^2 - dU. This is checked by finite
     differences on all 16 sign patterns. It has a positive common
     factor and the potential with a minus sign. That is the structure
     of signature_of_the_saddle.py S1, so S2 selects the same four
     Lorentzian cases.  The q = 0 tensor term vanishes identically with
     no tuning (composite_metric.py C4). The pole and G are those of the
     full form (composite_pole.py K8). The traceless part of the energy
     momentum tensor vanishes, and its trace is the constant Lambda. A
     c-number has no matrix element between one particle states, so
     Weinberg-Witten does not apply.

  L6 (the stress tensor is pure trace). For L = 2|det h|^(1/4) - U the
     canonical tensor is T^mu_nu = (U - |det h|^(1/4)) delta^mu_nu, from
     dL/dh = (1/2)|det h|^(1/4) h^-1. Its traceless part vanishes
     identically. Conservation makes U - |det h|^(1/4) a constant on
     shell, the integration constant Lambda.  A c-number has no matrix
     element between one particle states, so Weinberg-Witten does not
     apply.  The full form gives T = 0 identically (composite_metric.py
     C1).  Checked by finite differences.

RESULT. Given the axiom, the lane's three statements and locality, the
vacuum term has one closing form: the traceless form, with Lambda an
integration constant.  Its fixed density is case (iii) of the
classification, forced by the axiom.  A fixed density is held
structure, so the axiom and "nothing holds the field" cannot both hold
strictly, and the one density the axiom admits is this one.  The scope stated here (algebraic in U,
first order fields, D = 4, the gauge case by theta bar = 0) is reduced
in lambda_scope.py to premises the lane already holds:
  - the null Lagrangian theorem covers every local order;
  - translation invariance of the empty field vacuum removes every
    divergence built from the fields;
  - the gauge case is closed by the unnamed common mode, with
    theta bar = 0 as an independent backup;
  - global normalizations are fixed by the anchor;
  - every other nonlocal term is excluded by self-presence, the
    axiom's own locality.

Usage: python3 tests/probes/lambda_lane.py     (a second)
"""

import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import signature_of_the_saddle as sig  # noqa: E402

EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


def full(h, U):
    return math.sqrt(abs(np.linalg.det(h))) / U


def trace_free(h, U):
    return 2 * abs(np.linalg.det(h))**0.25 - U


def l1(rng, n=20, U=1.3, d=1e-4):
    worst_full, worst_tf = 0.0, 0.0
    for _ in range(n):
        g = rng.normal(size=(6, 4))
        h = g.T @ g
        d2 = lambda f: (f(h, U + d) - 2 * f(h, U) + f(h, U - d)) / d**2
        worst_full = max(worst_full, abs(d2(full) / (2 * full(h, U) / U**2) - 1))
        d1 = (trace_free(h, U + d) - trace_free(h, U - d)) / (2 * d)
        worst_tf = max(worst_tf, abs(d1 + 1) + abs(d2(trace_free)))
    return dict(full_second=worst_full, trace_free=worst_tf)


def random_closed_form(rng, n=6):
    """A constant coefficient 4-form on R^n (closed). Fully antisymmetric."""
    W = np.zeros((n,) * 4)
    for idx in itertools.combinations(range(n), 4):
        v = rng.normal()
        for p in itertools.permutations(range(4)):
            W[tuple(idx[i] for i in p)] = np.linalg.det(np.eye(4)[list(p)]) * v
    return W


def pullback(W, grad):
    """mu = eps^{mu nu rho sigma} W_ijkl d_mu phi^i d_nu phi^j d_rho phi^k d_sigma phi^l."""
    return float(np.einsum("abcd,ijkl,ia,jb,kc,ld->", EPS4, W, grad, grad, grad, grad))


def l2(rng, n=500):
    W = random_closed_form(rng)
    min_eig, rank3 = np.inf, 0.0
    for _ in range(n):
        grad = rng.normal(size=(6, 4))
        if abs(pullback(W, grad)) > 1e-6:
            min_eig = min(min_eig, np.linalg.eigvalsh(grad.T @ grad).min())
    for _ in range(100):
        grad = rng.normal(size=(6, 3)) @ rng.normal(size=(3, 4))
        rank3 = max(rank3, abs(pullback(W, grad)))
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    S = rng.normal(size=(6, 6))
    S = S + S.T
    A = rng.normal(size=(6, 6))
    A = A - A.T
    B = rng.normal(size=(4, 4))
    B = B - B.T

    def wick(G):
        # G[i, mu, j, nu] = <d_mu phi^i d_nu phi^j>; three pairings of four factors
        t = (np.einsum("iajb,kcld->ijklabcd", G, G)
             + np.einsum("iakc,jbld->ijklabcd", G, G)
             + np.einsum("iald,jbkc->ijklabcd", G, G))
        return float(np.einsum("abcd,ijkl,ijklabcd->", EPS4, W, t))
    lorentz = wick(np.einsum("ij,ab->iajb", S, eta))
    framed = wick(np.einsum("ij,ab->iajb", S, eta) + np.einsum("ij,ab->iajb", A, B))
    return dict(min_eig=float(min_eig), rank3=rank3, lorentz=abs(lorentz), framed=abs(framed))


def l5(rng, U=1.3, eps=1e-6):
    worst_kin = worst_pot = 0.0
    for pat in sig.PATTERNS:
        hv = np.diag(pat).astype(float)
        grad = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
        s = sig.gram(grad)
        fd = (trace_free(hv + eps * s, U) - trace_free(hv - eps * s, U)) / (2 * eps)
        kin = 0.5 * abs(np.linalg.det(hv))**0.25 * 2 * sum(
            pat[m] * np.sum(abs(grad[:, m])**2) for m in range(4))
        worst_kin = max(worst_kin, abs(fd / kin - 1))
        dU = 0.37
        fdU = (trace_free(hv, U + eps * dU) - trace_free(hv, U - eps * dU)) / (2 * eps)
        worst_pot = max(worst_pot, abs(fdU / (-dU) - 1))
    passing, _ = sig.hamiltonian_count()
    return dict(kin=worst_kin, pot=worst_pot, passing=len(passing))


def _matchings(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for k in range(len(rest)):
        for m in _matchings(rest[:k] + rest[k + 1:]):
            yield [(first, rest[k])] + m


def wick_field_dependent(C, P, M, D):
    """<C_ijklmn phi^m phi^n eps^abcd d_a phi^i d_b phi^j d_c phi^k d_d phi^l>, exact.

    P[m, n] = <phi phi>, M[m, i, a] = <phi^m d_a phi^i>, D[i, j, a, b] = <d_a phi^i d_b phi^j>.
    """
    slots = [("f", "m", None), ("f", "n", None), ("d", "i", "a"),
             ("d", "j", "b"), ("d", "k", "c"), ("d", "l", "e")]
    total = 0.0
    for match in _matchings(list(range(6))):
        ops, subs = [C, EPS4], ["ijklmn", "abce"]
        for p, q in match:
            sp, sq = slots[p], slots[q]
            if sp[0] == "f" and sq[0] == "f":
                ops.append(P); subs.append(sp[1] + sq[1])
            elif sp[0] == "f":
                ops.append(M); subs.append(sp[1] + sq[1] + sq[2])
            elif sq[0] == "f":
                ops.append(M); subs.append(sq[1] + sp[1] + sp[2])
            else:
                ops.append(D); subs.append(sp[1] + sq[1] + sp[2] + sq[2])
        total += float(np.einsum(",".join(subs) + "->", *ops, optimize=True))
    return total


def l2b(rng):
    n = 6
    W = random_closed_form(rng)
    Q = rng.normal(size=(n, n))
    C = np.einsum("ijkl,mn->ijklmn", W, Q + Q.T)      # coefficients quadratic in the fields
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    S = rng.normal(size=(n, n)); S = S + S.T
    G = rng.normal(size=(n, n)); G = G @ G.T
    lorentz = wick_field_dependent(C, G, np.zeros((n, n, 4)), np.einsum("ij,ab->ijab", S, eta))
    A = rng.normal(size=(n, n)); A = A - A.T
    B = rng.normal(size=(4, 4)); B = B - B.T
    v = rng.normal(size=4)
    framed = wick_field_dependent(C, G, np.einsum("mi,a->mia", rng.normal(size=(n, n)), v),
                                  np.einsum("ij,ab->ijab", S, eta) + np.einsum("ij,ab->ijab", A, B))
    return dict(lorentz=abs(lorentz), framed=abs(framed))


def l6(rng, U=1.3, d=1e-6):
    """Canonical T^mu_nu = sum_r dL/dG_{r mu} G_{r nu} - delta L for h = 2 G^T G."""
    worst_tf, worst_full = 0.0, 0.0
    for _ in range(10):
        G = rng.normal(size=(6, 4))
        for form, want in ((trace_free, lambda h: (U - abs(np.linalg.det(h))**0.25) * np.eye(4)),
                           (full, lambda h: np.zeros((4, 4)))):
            L = lambda g: form(2 * g.T @ g, U)
            dL = np.zeros_like(G)
            for r in range(6):
                for mu in range(4):
                    e = np.zeros_like(G); e[r, mu] = d
                    dL[r, mu] = (L(G + e) - L(G - e)) / (2 * d)
            T = dL.T @ G - np.eye(4) * L(G)
            dev = float(np.max(abs(T - want(2 * G.T @ G))))
            if form is trace_free:
                worst_tf = max(worst_tf, dev)
            else:
                worst_full = max(worst_full, dev)
    return dict(trace_free=worst_tf, full=worst_full)


def run(report=print, seed=7):
    rng = np.random.default_rng(seed)
    out = dict(l1=l1(rng), l2=l2(rng), l5=l5(rng), l2b=l2b(rng), l6=l6(rng))
    a, b, c = out["l1"], out["l2"], out["l5"]
    assert a["full_second"] < 1e-5 and a["trace_free"] < 1e-6
    assert b["min_eig"] > 0 and b["rank3"] < 1e-9
    assert b["lorentz"] < 1e-9 and b["framed"] > 1e-3
    assert c["kin"] < 1e-6 and c["pot"] < 1e-6 and c["passing"] == 4
    report(f"L1 full form: d2L/dU2 = 2 sqrt|det h|/U^3 to {a['full_second']:.1e}, so the constant "
           f"weighs histories. Traceless: dL/dU = -1, d2L/dU2 = 0 to {a['trace_free']:.1e}.")
    report(f"L2 a random closed 4-form on the six components: where its pullback is nonzero the "
           f"smallest eigenvalue of h is {b['min_eig']:.2e} > 0 (Euclidean). Gradients in three "
           f"directions give {b['rank3']:.1e}; Wick <mu> in a Lorentz invariant Gaussian state "
           f"{b['lorentz']:.1e}, with a frame fixing part {b['framed']:.2f}.")
    report("L3 no potential: no mass term, no self-bound sign, no lump.  L4 a fixed density "
           "with a quadratic kinetic term leaves c |det h|^(1/4): the traceless form.")
    report(f"L5 traceless first order term against (1/2)|det h_vac|^(1/4) sum |d phi|^2/sigma "
           f"to {c['kin']:.1e}, potential with a minus sign to {c['pot']:.1e}: the saddle "
           f"structure of the full form, {c['passing']} Lorentzian cases pass.")
    e, t = out["l2b"], out["l6"]
    assert e["lorentz"] < 1e-8 * max(1.0, e["framed"]) and e["framed"] > 1e-3
    assert t["trace_free"] < 1e-6 and t["full"] < 1e-6
    report(f"L2b coefficients quadratic in the fields, exact Wick over 15 pairings: Lorentz invariant "
           f"state {e['lorentz']:.1e}, frame fixing state {e['framed']:.3e}.")
    report(f"L6 canonical T of the traceless form against (U - |det h|^(1/4)) delta to "
           f"{t['trace_free']:.1e}. Full form T = 0 to {t['full']:.1e}.")
    report("RESULT: under the axiom the vacuum term has one closing form, the traceless form. Its "
           "fixed density is forced by the classification, and Lambda is an integration constant.")
    return out


if __name__ == "__main__":
    run()
