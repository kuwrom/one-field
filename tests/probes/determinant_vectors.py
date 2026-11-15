"""The carrier with a gauge field, written without a background metric.

THE QUESTION. composite_metric.py writes the carrier as the determinant
of its gradient bilinear, sqrt|det h| / U, with the metric g = h / U
eliminated. A gauge field needs the inverse metric to contract its field
strength, so the question is whether a closed form exists with F inside
the same determinant, and whether that form is the Maxwell theory on the
composite metric at quadratic order.

PREMISES. P1 the conventions of composite_metric.py C1: the auxiliary
form L = sqrt|g| (g^{mu nu} h_mu nu / 2 - U - F_mu nu F^mu nu / 4), with
F the field strength of an elementary vector and U > 0. P2 D = 4. P3 for
V3 only, F small against U, the expansion the whole construction uses
(arXiv:1610.08521 Sec. III).

THE FORM. For a symmetric h and an antisymmetric F in four dimensions
    det(h + beta F) = det h [1 + (beta^2/2) F_h^2] + beta^4 Pf(F)^2,
    F_h^2 = h^{mu rho} h^{nu sigma} F_mu nu F_rho sigma,
Pf is the Pfaffian, with Pf(F)^2 = (1/16) (F Ftilde)^2 |det h| and the
Lorentzian sign of arXiv:hep-th/9908105 Eq. 1.3 when Ftilde carries
1/sqrt|h|. Expanding the square root, a real beta gives the F^2 term the
sign of the 1/U term, which is the opposite of Maxwell's. Setting
beta^2 = -U fixes the sign, so the closed form is the determinant of a
Hermitian matrix,
    S = int sqrt|det(h_mu nu + i sqrt(U) F_mu nu)| / U.
It is real because the determinant of a Hermitian matrix is real. It
has the Dirac-Born-Infeld structure of the string action
(hep-th/9908105 Eq. 1.1, where the overall minus sign does the same
work). No overall sign
and no |U| is needed for U > 0. For several vectors or a nonabelian
field the determinant needs a trace prescription (the symmetrised trace,
hep-th/9908105 Sec. 5). Every prescription gives Yang-Mills at quadratic
order, so nothing below depends on it.

STATEMENTS (the numbers below check each one).

  V1 (identity). The four dimensional determinant identity above holds
     to machine precision for Euclidean and Lorentzian h, for real and
     imaginary beta, Pfaffian term included.
  V2 (Maxwell). To second order in F the Hermitian form is
     sqrt|det h| / U - (1/4) sqrt|det h| F_h^2: the Maxwell density on
     the composite metric, and U cancels from the F^2 term. The vector
     couples only to the conformal class of h (Weyl invariance of
     Maxwell in D = 4).
  V3 (elimination). Solving the algebraic metric equation of P1 exactly
     and evaluating the action on shell agrees with the Hermitian form at
     order F^2 and differs at order F^4 by
         sqrt|det h| U [-(1/32) (F_h^2)^2 + (1/2) Pf(F)^2 / det h],
     a contact term of the vector. At quadratic order both are metric
     free forms of one theory. The determinant form is closed. The
     eliminated form is a series in F^2 / U.
  V4 (invariances). The density is unchanged under A -> A + d chi and
     scales by |det M| under a general linear change of coordinates,
     h -> M^T h M, F -> M^T F M, U a scalar.
  V5 (the row). On a background h = c eta the quadratic form in A is the
     flat Maxwell form, independent of c and of U. The vector's one loop
     stiffness on the composite background is therefore the row already
     computed on a background metric in ledger_kernel.py (the VECTOR row,
     with its ghosts, F0/F2 -> -2). Nothing is recomputed here.

Usage: python3 tests/probes/determinant_vectors.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.optimize import root

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401

ETA = np.diag([1.0, -1.0, -1.0, -1.0])
IU = np.triu_indices(4)


def random_h(rng, lorentzian):
    A = rng.normal(size=(4, 4))
    return A.T @ (ETA if lorentzian else np.eye(4)) @ A


def random_f(rng):
    B = rng.normal(size=(4, 4))
    return B - B.T


def pfaffian(F):
    return F[0, 1] * F[2, 3] - F[0, 2] * F[1, 3] + F[0, 3] * F[1, 2]


def f_sq(h, F):
    """F_h^2 = h^{mu rho} h^{nu sigma} F_mu nu F_rho sigma."""
    hi = np.linalg.inv(h)
    return float(np.trace(hi @ F @ hi @ F.T))


def identity_residual(h, F, beta):
    """V1: det(h + beta F) against the closed form, relative."""
    lhs = np.linalg.det(h + beta * F)
    dh = np.linalg.det(h)
    rhs = dh * (1 + 0.5 * beta**2 * f_sq(h, F)) + beta**4 * pfaffian(F)**2
    return abs(lhs - rhs) / (abs(dh) + abs(lhs))


def hermitian_form(h, F, U):
    """sqrt|det(h + i sqrt U F)| / U."""
    return math.sqrt(abs(np.linalg.det(h + 1j * math.sqrt(U) * F).real)) / U


def maxwell(h, F):
    """-(1/4) sqrt|det h| F_h^2, the Maxwell density on g = h / U."""
    return -0.25 * math.sqrt(abs(np.linalg.det(h))) * f_sq(h, F)


def lagrangian(g, h, U, F):
    """P1: sqrt|g| (g.h/2 - U - F.F/4) with indices raised by g."""
    gi = np.linalg.inv(g)
    return math.sqrt(abs(np.linalg.det(g))) * (
        0.5 * np.trace(gi @ h) - U - 0.25 * np.trace(gi @ F @ gi @ F.T))


def metric_equation(g, h, U, F):
    """Stationarity in g: h - F g^{-1} F^T - g (g.h/2 - U - F.F/4) = 0."""
    gi = np.linalg.inv(g)
    ff = np.trace(gi @ F @ gi @ F.T)
    return h - F @ gi @ F.T - g * (0.5 * np.trace(gi @ h) - U - 0.25 * ff)


def eliminate(h, U, F):
    """V3: solve the metric equation from g = h/U and return the on shell
    action and the residual."""
    def unpack(x):
        g = np.zeros((4, 4))
        g[IU] = x
        return g + g.T - np.diag(np.diag(g))

    def residual(x):
        return metric_equation(unpack(x), h, U, F)[IU]
    sol = root(residual, (h / U)[IU], method="hybr", tol=1e-14)
    g = unpack(sol.x)
    return lagrangian(g, h, U, F), float(np.max(abs(metric_equation(g, h, U, F))))


def contact_fit(rng, U, eps=0.05, samples=12):
    """V3: fit (L_exact - L_det) / (sqrt|det h| U eps^4) to
    c1 (F_h^2)^2 + c2 Pf(F)^2 / det h over random Lorentzian h and F."""
    rows, rhs, resid = [], [], 0.0
    for _ in range(samples):
        h, F = random_h(rng, True), random_f(rng)
        F *= 1.0 / math.sqrt(abs(f_sq(h, F)))
        d = abs(np.linalg.det(h))**0.5
        ex, r = eliminate(h, U, eps * F)
        resid = max(resid, r)
        diff = ex - hermitian_form(h, eps * F, U)
        rows.append([f_sq(h, F)**2, pfaffian(F)**2 / np.linalg.det(h)])
        rhs.append(diff / (d * U * eps**4))
    rows, rhs = np.array(rows), np.array(rhs)
    coef, *_ = np.linalg.lstsq(rows, rhs, rcond=None)
    misfit = float(np.max(abs(rows @ coef - rhs)) / np.max(abs(rhs)))
    return coef, misfit, resid


def quadratic_in_f(h, F, U, eps=1e-3):
    """Second order coefficient of the Hermitian form in F by central
    differences, with the fourth order removed by Richardson."""
    def d2(e):
        return (hermitian_form(h, e * F, U) - 2 * hermitian_form(h, 0 * F, U)
                + hermitian_form(h, -e * F, U)) / e**2
    return 0.5 * (4 * d2(eps) - d2(2 * eps)) / 3


def run(report=print):
    rng = np.random.default_rng(7)
    out = {}

    # V1 the determinant identity
    worst = 0.0
    for lor in (False, True):
        for beta in (0.7, 1j * 0.9, 1.3j):
            for _ in range(5):
                worst = max(worst, identity_residual(random_h(rng, lor), random_f(rng), beta))
    out["identity"] = worst
    report(f"V1 det(h + beta F) = det h [1 + (beta^2/2) F_h^2] + beta^4 Pf(F)^2: "
           f"worst relative residual {worst:.1e} (Euclidean and Lorentzian h, "
           f"real and imaginary beta).")

    # V2 the Maxwell term and its U independence
    err_max, err_u = 0.0, 0.0
    for lor in (False, True):
        for _ in range(4):
            h, F = random_h(rng, lor), random_f(rng)
            q1, q2 = quadratic_in_f(h, F, 0.8), quadratic_in_f(h, F, 3.1)
            err_max = max(err_max, abs(q1 / maxwell(h, F) - 1))
            err_u = max(err_u, abs(q1 / q2 - 1))
    out["maxwell"], out["u_free"] = err_max, err_u
    report(f"V2 second order of sqrt|det(h + i sqrt U F)|/U against -(1/4) sqrt|det h| "
           f"F_h^2: relative error {err_max:.1e}. U = 0.8 against U = 3.1: {err_u:.1e}.")

    # V3 exact elimination against the closed form
    h, F, U = random_h(rng, True), random_f(rng), 1.4
    F *= 1.0 / math.sqrt(abs(f_sq(h, F)))
    d = math.sqrt(abs(np.linalg.det(h)))
    order2 = []
    for eps in (0.2, 0.1, 0.05):
        ex, _ = eliminate(h, U, eps * F)
        order2.append((ex - hermitian_form(h, eps * F, U)) / (d * eps**2))
    coef, misfit, resid = contact_fit(rng, U)
    out["order2"], out["contact"], out["contact_misfit"], out["residual"] = order2, coef, misfit, resid
    report(f"V3 on shell action minus the Hermitian form, over sqrt|det h| eps^2 at "
           f"eps = 0.2, 0.1, 0.05: {order2[0]:.2e}, {order2[1]:.2e}, {order2[2]:.2e} "
           f"(vanishing as eps^2: agreement at order F^2). The F^4 difference is "
           f"sqrt|det h| U [{coef[0]:+.5f} (F_h^2)^2 {coef[1]:+.5f} Pf^2/det h], "
           f"fit misfit {misfit:.1e}, metric equation residual {resid:.1e} "
           f"(-1/32 = {-1 / 32:.5f}, 1/2).")

    # V4 gauge and general linear invariance
    B, Q, x = rng.normal(size=(4, 4)), rng.normal(size=(4, 4)), rng.normal(size=4)
    Q = Q + Q.T                                    # chi = x Q x / 2, d chi = Q x
    A1 = lambda y: B @ y                           # noqa: E731
    A2 = lambda y: B @ y + Q @ y                   # noqa: E731
    F1 = B.T - B                                   # F_mu nu = d_mu A_nu - d_nu A_mu
    F2 = (B + Q).T - (B + Q)
    h, U = random_h(rng, True), 2.0
    gauge = abs(hermitian_form(h, F1, U) - hermitian_form(h, F2, U))
    M = rng.normal(size=(4, 4))
    cov = abs(hermitian_form(M.T @ h @ M, M.T @ F1 @ M, U)
              / (abs(np.linalg.det(M)) * hermitian_form(h, F1, U)) - 1)
    out["gauge"], out["covariance"] = gauge, cov
    report(f"V4 A -> A + d chi: change {gauge:.1e} (A2 - A1 at a point: "
           f"{np.max(abs(A2(x) - A1(x) - Q @ x)):.1e}). General linear map: density "
           f"scales by |det M| to {cov:.1e}.")

    # V5 flat background: the quadratic form is the flat Maxwell form for every c
    worst = 0.0
    F = random_f(rng)
    flat = maxwell(ETA, F)
    for c in (0.3, 1.0, 4.0):
        for U in (0.5, 2.0):
            worst = max(worst, abs(quadratic_in_f(c * ETA, F, U) / flat - 1))
    out["flat"] = worst
    report(f"V5 on h = c eta the second order form in A equals the flat Maxwell form "
           f"for c = 0.3, 1, 4 and U = 0.5, 2 to {worst:.1e}: the vector row of "
           f"ledger_kernel.py carries over unchanged.")
    report("RESULT: a closed metric free form with a gauge field exists, the "
           "determinant of h + i sqrt(U) F over U. At quadratic order it is "
           "Maxwell on the composite metric, so the vector row of the ledger is "
           "the heat trace row. The forms differ by a contact term at order F^4.")
    return out


if __name__ == "__main__":
    run()
