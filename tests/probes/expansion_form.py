"""The expansion's form: the traceless equations on a homogeneous isotropic metric.

THE QUESTION.  The lane's field equations are the traceless ones,
"with Lambda an integration constant" (lambda_lane.py), and the amount
is read as "a count, conserved per comoving volume, whose weight per
count does not depend on the volume" (dark_sector.py, THE DUST
READING).  What do the two give for the expansion of a homogeneous
isotropic universe, and what do they leave to the state?

PREMISES.  P1 the traceless equations
    R_mn - (1/4) g_mn R = 8 pi G (T_mn - (1/4) g_mn T)
with the G of the static and radiated fields (post_newtonian.py).
P2 the metric -dt^2 + a(t)^2 (dr^2/(1 - k r^2) + r^2 dOmega^2) and a
perfect fluid source of density rho and pressure p.  P3 the source is
conserved, d rho/dt = -3 H (rho + p). For the amount this is the
conserved count of the dust reading. For the nameable rows it is their
own conservation.  The
traceless equations do not imply P3 as the Einstein equations do.
Their divergence gives d(R + 8 pi G T) = 0 only when P3 holds.

STATEMENTS (the numbers below check each one).

  E1 (one equation).  Every component of P1 on P2 is the single
     equation dH/dt - k/a^2 = -4 pi G (rho + p).  None of them is the
     rate equation. The traceless equations fix dH/dt and leave the
     constant in the rate equation free.
  E2 (the constant).  With P3 the combination
     3 H^2 + 3 k/a^2 - 8 pi G rho has vanishing time derivative on E1.
     It is a constant of the motion, Lambda, and
     R + 8 pi G T = 4 Lambda.  The rate equation
     H^2 + k/a^2 = (8 pi G/3) rho + Lambda/3 is the first integral.
  E3 (the constant is initial data).  Integrated numerically from two
     initial rates with the same content, E1 keeps two different
     values of Lambda, each fixed by its initial data to 1e-9. The
     equations have no preferred value.
  E4 (its equation of state).  Treated as a density, the constant
     does not dilute: w = -1 exactly, at every redshift.  The count
     dilutes as a^-3, the nameable matter as a^-3 and radiation as
     a^-4, so the rate is H^2 = H0^2 (Or a^-4 + Om a^-3 + Ok a^-2 + OL).
     The ratio of the amount to the baryons is the same at every epoch,
     the 2 pi - 1 of dark_sector.py.
  E5 (what fails).  A source that is not conserved, d rho/dt =
     -3 H (rho + p) + Q, gives dLambda/dt = -8 pi G Q. There is then no
     constant, and the vacuum term moves.  An amount with pressure
     p = w rho dilutes as a^(-3(1 + w)), and its ratio to the baryons
     drifts by (1 + z)^(3 w). At the CMB's own bound w = 0.0024
     (arXiv:1601.05097) the ratio moves by 5.2 percent between
     z = 1100 and today.

RESULT.  The traceless equations and the conserved count give the
expansion the form of a constant vacuum term with pressureless matter.
That is the rate equation with Lambda as its integration constant,
w = -1, and a dark to baryon ratio that does not move.  The form does
not fix Lambda's value or the spatial curvature k. Those are data of
the state (notes/open_derivation_tasks.md T3), and the rate today is
an anchor (root.py).  A measured w(z) away from -1, or a ratio that
drifts, falsifies it.

Scope: homogeneous isotropic metric, perfect fluids, the background
expansion. Perturbations are treated in growth_form.py.

Usage: python3 tests/probes/expansion_form.py     (a few seconds)
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

W_BOUND = 0.0024                         # CMB upper bound on the dark component's w (arXiv:1601.05097)
Z_CMB = 1100.0


def trace_free_components():
    """Components of P1 on P2, each divided by its source factor: (name, left - right)."""
    t, rr, th, ph, k, G = sp.symbols("t r theta phi k G", real=True)
    a = sp.Function("a", positive=True)(t)
    rho, p = sp.Function("rho")(t), sp.Function("p")(t)
    x = (t, rr, th, ph)
    g = sp.diag(-1, a**2 / (1 - k * rr**2), a**2 * rr**2, a**2 * rr**2 * sp.sin(th)**2)
    ginv = g.inv()
    Gam = [[[sum(ginv[i, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d]))
                 for d in range(4)) / 2 for c in range(4)] for b in range(4)] for i in range(4)]

    def ricci(b, c):
        out = 0
        for i in range(4):
            out += sp.diff(Gam[i][b][c], x[i]) - sp.diff(Gam[i][b][i], x[c])
            out += sum(Gam[i][i][d] * Gam[d][b][c] - Gam[i][c][d] * Gam[d][b][i] for d in range(4))
        return sp.simplify(out)

    Ric = sp.Matrix(4, 4, lambda b, c: ricci(b, c))
    R = sp.simplify(sum(ginv[b, c] * Ric[b, c] for b in range(4) for c in range(4)))
    u = sp.Matrix([-1, 0, 0, 0])                                  # u_mu of the comoving fluid
    T = sp.Matrix(4, 4, lambda b, c: (rho + p) * u[b] * u[c] + p * g[b, c])
    trT = sp.simplify(sum(ginv[b, c] * T[b, c] for b in range(4) for c in range(4)))
    E = sp.Matrix(4, 4, lambda b, c: sp.simplify(
        Ric[b, c] - g[b, c] * R / 4 - 8 * sp.pi * G * (T[b, c] - g[b, c] * trT / 4)))
    H = sp.diff(a, t) / a
    one = sp.diff(H, t) - k / a**2 + 4 * sp.pi * G * (rho + p)    # E1
    return dict(E=E, g=g, R=R, trT=trT, one=one, H=H, symbols=(t, a, rho, p, k, G))


def constant_of_motion(parts):
    """d/dt of 3 H^2 + 3 k/a^2 - 8 pi G rho on E1, with the source's conservation broken by Q."""
    t, a, rho, p, k, G = parts["symbols"]
    H = parts["H"]
    Q = sp.Symbol("Q")
    lam = 3 * H**2 + 3 * k / a**2 - 8 * sp.pi * G * rho
    addot = sp.solve(parts["one"], sp.diff(a, t, 2))[0]
    rate = sp.diff(lam, t).subs(sp.diff(a, t, 2), addot).subs(sp.diff(rho, t), -3 * H * (rho + p) + Q)
    trace = (parts["R"] + 8 * sp.pi * G * parts["trT"]).subs(sp.diff(a, t, 2), addot) - 4 * lam
    return sp.simplify(rate), sp.simplify(trace), Q, G


def integrate(H0, k=0.0, matter=0.3, radiation=1e-4, a_end=3.0):
    """E1 integrated from a = 1 in units 8 pi G/3 = 1. Returns Lambda/3 along the solution."""
    rho = lambda a: matter / a**3 + radiation / a**4
    rho_plus_p = lambda a: matter / a**3 + (4.0 / 3.0) * radiation / a**4

    def rhs(_, y):
        a, H = y
        return [a * H, k / a**2 - 1.5 * rho_plus_p(a)]

    sol = solve_ivp(rhs, (0.0, 50.0), [1.0, H0], rtol=1e-11, atol=1e-13, dense_output=True,
                    events=lambda _, y: y[0] - a_end)
    a, H = sol.y
    lam = H**2 + k / a**2 - rho(a)
    return a, H, lam


def run(report=print):
    out = {}

    # E1: one equation.
    parts = trace_free_components()
    E, g, one = parts["E"], parts["g"], parts["one"]
    t, a, rho, p, k, G = parts["symbols"]
    ratios = []
    for i in range(4):
        for j in range(4):
            if i != j:
                assert E[i, j] == 0
        ratios.append(sp.simplify(E[i, i] / one))
    assert all(sp.simplify(sp.diff(q, sp.diff(a, t, 2))) == 0 and q != 0 for q in ratios)
    assert sp.simplify(ratios[0] + sp.Rational(3, 2)) == 0
    assert all(sp.simplify(ratios[i] + g[i, i] / 2) == 0 for i in (1, 2, 3))
    out.update(components=4, independent=1)
    report("E1 one equation: the four diagonal components of the traceless equations are "
           "dH/dt - k/a^2 + 4 pi G (rho + p) times -3/2 and -g_ii/2. The off diagonal ones vanish. "
           "None is the rate equation")

    # E2, E5a: the constant, and its rate when the source is not conserved.
    rate, trace, Q, Gs = constant_of_motion(parts)
    assert sp.simplify(rate.subs(Q, 0)) == 0 and trace == 0
    assert sp.simplify(rate + 8 * sp.pi * Gs * Q) == 0
    out.update(constant=True, rate_with_Q="-8 pi G Q")
    report("E2 constant: d/dt (3 H^2 + 3 k/a^2 - 8 pi G rho) = 0 on E1 with a conserved source, and "
           "R + 8 pi G T is four times it")

    # E3: the constant is initial data.
    runs = {H0: integrate(H0) for H0 in (1.0, 1.2)}
    lams = {H0: (lam[0], float(np.max(np.abs(lam - lam[0])))) for H0, (_, _, lam) in runs.items()}
    assert all(drift < 1e-9 for _, drift in lams.values())
    assert abs(lams[1.0][0] - (1.0 - 0.3 - 1e-4)) < 1e-12 and lams[1.2][0] - lams[1.0][0] > 0.4
    curved = integrate(1.0, k=0.2)
    assert float(np.max(np.abs(curved[2] - curved[2][0]))) < 1e-9
    out.update(lambdas={H0: v[0] for H0, v in lams.items()}, drift=max(v[1] for v in lams.values()))
    report("E3 initial data: from H0 = 1.0 and 1.2 with the same content, Lambda/3 = "
           + " and ".join(f"{v[0]:.4f}" for v in lams.values())
           + f", each constant to {max(v[1] for v in lams.values()):.0e} along the solution. With k = 0.2 "
           "it is constant as well")

    # E4: w = -1 and the rate's form.
    a_run, H_run, lam_run = runs[1.0]
    form = 1e-4 / a_run**4 + 0.3 / a_run**3 + lam_run[0]
    w_eff = -1.0 - np.gradient(np.log(np.abs(lam_run)), np.log(a_run)) / 3.0
    assert float(np.max(np.abs(H_run**2 / form - 1.0))) < 1e-9 and float(np.max(np.abs(w_eff + 1.0))) < 1e-7
    out.update(form_error=float(np.max(np.abs(H_run**2 / form - 1.0))), w=-1.0,
               w_error=float(np.max(np.abs(w_eff + 1.0))))
    report(f"E4 equation of state: H^2 is Or a^-4 + Om a^-3 + OL to {out['form_error']:.0e} from a = 1 to 3, "
           f"and the constant's w is -1 to {out['w_error']:.0e}. The amount and the baryons both dilute as a^-3")

    # E5b: an amount with pressure.
    drift = (1.0 + Z_CMB)**(3 * W_BOUND) - 1.0
    assert 0.05 < drift < 0.054 and (1.0 + Z_CMB)**(3 * 0.0) - 1.0 == 0.0
    out.update(ratio_drift_at_bound=drift)
    report(f"E5 what fails: a source with exchange Q gives dLambda/dt = -8 pi G Q. An amount with "
           f"w = {W_BOUND} moves its ratio to the baryons by {100 * drift:.1f} percent between z = {Z_CMB:.0f} "
           "and today")

    report("RESULT: the traceless equations with a conserved count give the rate equation with Lambda its "
           "integration constant, w = -1, and a dark to baryon ratio that does not move. Lambda's value and k "
           "are data of the state.")
    return out


if __name__ == "__main__":
    run()
