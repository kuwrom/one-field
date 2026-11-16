"""The post-Newtonian gamma and beta, and one G for static and radiated fields.

THE QUESTION. The induced kernel gives Newton's law, and its spin
content gives gamma (composite_pole.py K5). The two derivative term of
the induced action is the Einstein-Hilbert density beyond quadratic
order, for every row (induced_cubic.py, induced_cubic_rows.py). What
are gamma and beta for that action, as numbers, and does a radiated
field couple with the G of the static one?

PREMISES. P1 the action int sqrt g R through cubic order, in vacuum
outside a static spherical source. P2 the isotropic post-Newtonian
form, U = G M / r:
    g_tt = -(1 - 2 U + 2 beta U^2),
    g_ij = (1 + 2 gamma U + (3/2) delta U^2) delta_ij.
P3 the kernel of composite_pole.py under the proper time cutoff, at
spacelike and at timelike momenta.

STATEMENTS (the numbers below check each one).

  N1 (first order). The vacuum equations at first order in M hold for
     gamma = 1 alone: the kernel's F0/F2 = -2 is the same statement,
     gamma = (1 - F2/F0)/(2 + F2/F0). F0/F2 = -2.2 would give
     gamma = 0.941, a sensitivity check and not a rival reading: the
     kernel's q^2 term is local only at -2 (composite_pole.py K9).
  N2 (second order). At second order in M, with gamma = 1, the
     equations hold for beta = 1 and delta = 1 alone, the expansion of
     the Schwarzschild metric. beta depends on the cubic term of the
     action. With the Einstein cubic term the induced action gives
     beta = 1.
  N3 (radiation). The kernel depends on the momentum through q^2
     alone (saddle_stabilizer.py Z1). It is one analytic function on
     both sides of the light cone. At |q^2| = 1e-4 m^2, 4 F2/q^2 over the
     stiffness is 1 - 1.5e-5 at spacelike and 1 + 1.5e-5 at timelike
     momentum, the same slope with the sign of q^2, and their mean is
     the stiffness to 1e-8. The pole's residue is a single number, so
     a radiated wave couples with the G of Newton's law.

RESULT. For the induced action gamma = beta = 1, and the spin-2 pole
that carries Newton's law carries radiation with the same G.

Scope: vacuum outside a static spherical source to second order.
Radiation at the level of the pole's residue.

Usage: python3 tests/probes/post_newtonian.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import composite_pole as cp  # noqa: E402


def ricci_series(order=2):
    """Ricci components of the P2 metric as polynomials in M, to the given order."""
    t, r, th, ph, M = sp.symbols("t r theta phi M", positive=True)
    beta, gamma, delta = sp.symbols("beta gamma delta")
    U = M / r
    A = 1 - 2 * U + 2 * beta * U**2
    Bf = 1 + 2 * gamma * U + sp.Rational(3, 2) * delta * U**2
    x = (t, r, th, ph)
    g = sp.diag(-A, Bf, Bf * r**2, Bf * r**2 * sp.sin(th)**2)
    ginv = g.inv()
    Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d]))
                 for d in range(4)) / 2 for c in range(4)] for b in range(4)] for a in range(4)]

    def ricci(b, c):
        out = 0
        for a in range(4):
            out += sp.diff(Gam[a][b][c], x[a]) - sp.diff(Gam[a][b][a], x[c])
            out += sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(4))
        return sp.series(sp.simplify(out), M, 0, order + 1).removeO()
    return {"tt": ricci(0, 0), "rr": ricci(1, 1), "thth": ricci(2, 2)}, (M, r, beta, gamma, delta)


def run(report=print):
    out = {}
    ric, (M, r, beta, gamma, delta) = ricci_series()

    # N1: first order in M.
    first = [sp.simplify(e.coeff(M, 1)) for e in ric.values()]
    sol1 = sp.solve([e for e in first if e != 0], [gamma], dict=True)
    assert sol1 == [{gamma: 1}]
    kernel_gamma = lambda f0_over_f2: (1 - 1 / f0_over_f2) / (2 + 1 / f0_over_f2)
    assert abs(kernel_gamma(-2.0) - 1) < 1e-12 and abs(kernel_gamma(-2.2) - 0.9412) < 1e-4
    out.update(gamma=1, gamma_moved=kernel_gamma(-2.2))
    report(f"N1 first order: the vacuum equations fix gamma = {sol1[0][gamma]}. The kernel's F0/F2 = -2 gives "
           f"{kernel_gamma(-2.0):.4f}, and -2.2 would give {kernel_gamma(-2.2):.4f}")

    # N2: second order in M at gamma = 1.
    second = [sp.simplify(e.coeff(M, 2).subs(gamma, 1)) for e in ric.values()]
    sol2 = sp.solve(second, [beta, delta], dict=True)
    assert sol2 == [{beta: 1, delta: 1}]
    off = [sp.simplify(e.subs({beta: sp.Rational(11, 10), delta: 1})) for e in second]
    assert any(e != 0 for e in off)
    out.update(beta=1, delta=1)
    report(f"N2 second order: beta = {sol2[0][beta]}, delta = {sol2[0][delta]}. Beta = 1.1 leaves a residual")

    # N3: one stiffness on both sides of the cone.
    J, _, stiffness = cp.proper_time(1.0, 1.0)
    s = 1e-4
    space = 4 * cp.form_factors(np.array([math.sqrt(s), 0.0, 0.0, 0.0]), J)[0] / s
    time_like = 4 * cp.timelike_ratio(J, s)[0]
    lo, hi = space / stiffness - 1, float(np.real(time_like)) / stiffness - 1
    assert -2e-5 < lo < -1e-5 and 1e-5 < hi < 2e-5 and abs(lo + hi) < 2e-8
    out.update(spacelike=lo, timelike=hi)
    report(f"N3 radiation: 4 F2/q^2 over the stiffness is 1 {lo:+.2e} at spacelike and 1 {hi:+.2e} at timelike "
           f"momentum. Their mean is the stiffness to {abs(lo + hi) / 2:.0e}")
    report("RESULT: gamma = beta = 1 for the induced action, and radiation couples with Newton's G.")
    return out


if __name__ == "__main__":
    run()
