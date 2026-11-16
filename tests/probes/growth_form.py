"""Linear growth on the expansion's form: what the count reading gives for structure.

THE QUESTION.  The traceless equations with a conserved count give
the rate equation with Lambda its integration constant
(expansion_form.py), and the amount is read as a count with "p = 0"
and "c_s^2 = 0" (dark_sector.py, THE DUST READING).  How does a small
excess of the count grow on that background, on which scales, and what
does the ratio 2 pi - 1 of the amount to the baryons fix?

PREMISES.  P1 expansion_form.py E1 and E2: for pressureless matter the
acceleration is (d^2a/dt^2)/a = -(4 pi G/3) rho + Lambda/3, with Lambda
one constant of the whole solution, R + 8 pi G T = 4 Lambda at every
point.  P2 a spherical region of the count evolves as its own solution
of the same equations: outside a spherical source the static vacuum
solution is unique (strong_field.py F2), so the region does not feel
what lies around it.  P3 the count is conserved, rho a^3 constant in
each region, and has no pressure and no sound speed.  P4 the amount is
2 pi - 1 times the baryons, a fraction f = 1 - 1/(2 pi) of the matter
(dark_sector.py).  P5 for the numbers: a spatially flat background
with matter fraction 0.315 today, a datum of the state.

STATEMENTS (the numbers below check each one).

  D1 (the growth equation).  A region with scale factor a (1 - s) holds
     the density rho (1 + 3 s).  P1 for the region minus P1 for the
     background, to first order in s, is
         d^2 delta/dt^2 + 2 H d delta/dt = 4 pi G rho delta,
     delta = 3 s.  Lambda has dropped out.  A single constant, the
     same in the region and around it, carries no perturbation.  It
     enters only through H.
  D2 (no scale).  The equation contains no wavenumber: with c_s = 0 every
     scale grows by the same factor D(a).  A sound speed adds
     c_s^2 k^2/a^2 to the restoring side and stops growth below the
     Jeans length 2 pi c_s/sqrt(4 pi G rho): at the CMB's bound
     c_s^2 = 3.2e-6 (arXiv:1601.05097) that is 73 Mpc today.
  D3 (the growth factor).  In the matter era D is proportional to a.
     On the form of expansion_form.py E4 the rate f = d ln D/d ln a
     today is the matter fraction to the power 0.55, between 0.553 and
     0.557 for fractions from 0.25 to 0.35: the growth index of a
     constant vacuum term.
  D4 (what 2 pi - 1 fixes).  While the baryons are held smooth the
     amount alone sources its growth, and D goes as a^p with
     p = (sqrt(1 + 24 f) - 1)/4 = 0.9006 for f = 1 - 1/(2 pi).  The
     other two readings of boundary_whole.py B4 give a dark ratio R
     infinite (the remainder alone) or zero (the local share alone),
     so f = R/(1 + R) is 1 or 0 and the exponent 1 or 0.  Released at
     a = 1/1100 with no excess of their own, the baryons reach 0.97 of
     the amount's excess by a = 1/10 and 0.996 today.  There is one
     growing mode, and the ratio 2 pi - 1 holds in the structures as in
     the mean.
  D5 (what fails).  A count that is not conserved has no growth
     equation of this form (expansion_form.py E5).  The vacuum term is
     the traceless form's, the one lane that closes (lambda_lane.py),
     and its Lambda is one constant of the whole solution, so it has
     no excess of its own (D1).

RESULT.  The count grows by one scale free factor, proportional to a
in the matter era and with growth index 0.55 under the constant vacuum
term.  Lambda carries no perturbation. The amount's share 1 - 1/(2 pi)
sets the exponent 0.9006 while the baryons are held, and the baryons
then join the same mode.

Scope: linear order, pressureless matter, the separate universe
construction, the growth of a given initial perturbation.

Usage: python3 tests/probes/growth_form.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boundary_whole as bw  # noqa: E402

OMEGA_M = 0.315                         # the matter fraction today used for the numbers (P5)
CS2_BOUND = 3.2e-6                      # CMB bound on the dark component's sound speed squared, c = 1
H0_OVER_C = 67.4 / 299792.458           # per Mpc, the engine's H0 anchor (root.py)
F_AMOUNT = 1.0 - 1.0 / (2 * math.pi)    # the amount's share of the matter
A_RELEASE = 1.0 / 1100.0


def growth_equation():
    """P1 for a region a (1 - s) minus the background, to first order: the residual and its Lambda dependence."""
    t, G, Lam, eps = sp.symbols("t G Lambda epsilon", positive=True)
    a = sp.Function("a", positive=True)(t)
    rho = sp.Function("rho", positive=True)(t)
    s = sp.Function("s")(t)
    acc = lambda scale, dens: sp.diff(scale, t, 2) / scale + 4 * sp.pi * G * dens / 3 - Lam / 3
    region = acc(a * (1 - eps * s), rho * (1 + 3 * eps * s))
    first = sp.diff(region, eps).subs(eps, 0)                       # first order in the excess
    H = sp.diff(a, t) / a
    background = sp.solve(acc(a, rho), sp.diff(a, t, 2))[0]
    first = sp.simplify(first.subs(sp.diff(a, t, 2), background))
    delta = 3 * s
    target = -(sp.diff(delta, t, 2) + 2 * H * sp.diff(delta, t) - 4 * sp.pi * G * rho * delta) / 3
    return sp.simplify(first - target), first.has(Lam)


def growth(omega_m, a0=1e-3, source=lambda a: 1.0, n=400):
    """D(a) and d ln D/d ln a on the flat form with a constant vacuum term, from the matter era mode at a0."""
    def rhs(x, y):
        a = math.exp(x)
        e2 = omega_m / a**3 + 1.0 - omega_m
        om = omega_m / a**3 / e2
        return [y[1], -(2.0 - 1.5 * om) * y[1] + 1.5 * om * source(a) * y[0]]
    xs = np.linspace(math.log(a0), 0.0, n)
    sol = solve_ivp(rhs, (xs[0], 0.0), [a0, a0], t_eval=xs, rtol=1e-10, atol=1e-14)
    return np.exp(xs), sol.y[0], sol.y[1] / sol.y[0]


def held_exponent(f):
    """Matter era growth exponent of a share f of the matter while the rest is held smooth."""
    return (math.sqrt(1.0 + 24.0 * f) - 1.0) / 4.0


def release(omega_m=OMEGA_M, f=F_AMOUNT, a0=A_RELEASE):
    """Amount and baryon excesses after the baryons are released smooth at a0. Returns a, delta_c, delta_b."""
    p = held_exponent(f)

    def rhs(x, y):
        a = math.exp(x)
        e2 = omega_m / a**3 + 1.0 - omega_m
        om = omega_m / a**3 / e2
        src = 1.5 * om * (f * y[0] + (1.0 - f) * y[2])
        drag = 2.0 - 1.5 * om
        return [y[1], -drag * y[1] + src, y[3], -drag * y[3] + src]
    xs = np.linspace(math.log(a0), 0.0, 600)
    sol = solve_ivp(rhs, (xs[0], 0.0), [1.0, p, 0.0, 0.0], t_eval=xs, rtol=1e-10, atol=1e-14)
    return np.exp(xs), sol.y[0], sol.y[2]


def _silent(*args, **kwargs):
    pass


def run(report=print):
    out = {}

    # D1: the growth equation, and Lambda's absence.
    residual, holds_lambda = growth_equation()
    assert residual == 0 and not holds_lambda
    out.update(growth_equation=True, lambda_perturbs=False)
    report("D1 growth equation: the region minus the background is delta'' + 2 H delta' - 4 pi G rho delta = 0 "
           "to first order; Lambda does not appear")

    # D2: no scale, and the Jeans length a sound speed would give.
    k_jeans = math.sqrt(1.5 * OMEGA_M) * H0_OVER_C / math.sqrt(CS2_BOUND)
    jeans = 2 * math.pi / k_jeans
    assert 70 < jeans < 76
    out.update(jeans_Mpc=jeans, jeans_at_zero_sound_speed=0.0)
    report(f"D2 no scale: the equation contains no wavenumber at c_s = 0; a sound speed at the CMB bound gives a "
           f"Jeans length of {jeans:.0f} Mpc today")

    # D3: the growth factor and its index.
    a, D, f = growth(1.0)
    eds = float(np.max(np.abs(D / a - 1.0)))
    index = {}
    for om in (0.25, OMEGA_M, 0.35):
        _, _, rate = growth(om)
        index[om] = math.log(rate[-1]) / math.log(om)
    _, D_L, f_L = growth(OMEGA_M)
    assert eds < 1e-8 and all(0.553 < g < 0.557 for g in index.values())
    out.update(matter_era_error=eds, growth_index=index, rate_today=float(f_L[-1]),
               suppression=float(D_L[-1]))
    report(f"D3 growth factor: D = a in the matter era to {eds:.0e}. With the constant vacuum term the rate today "
           f"is {f_L[-1]:.4f} at a matter fraction {OMEGA_M}, index "
           + ", ".join(f"{g:.4f}" for g in index.values()) + " for fractions 0.25, 0.315, 0.35")

    # D4: the amount's share, the held exponent, and the baryons' release.
    p = held_exponent(F_AMOUNT)
    a_h, D_h, f_h = growth(1.0, source=lambda _: F_AMOUNT)
    numeric = float(f_h[-1])
    ratios = bw.run(report=_silent)["dark_readings"]                  # boundary_whole.py B4's three readings
    shares = {k: 1.0 if math.isinf(R) else R / (1.0 + R) for k, R in ratios.items()}
    held = {k: held_exponent(f) for k, f in shares.items()}
    a_r, dc, db = release()
    at_tenth = float(np.interp(0.1, a_r, db / dc))
    today = float(db[-1] / dc[-1])
    assert abs(p - 0.9006) < 1e-4 and abs(numeric - p) < 1e-6
    assert abs(shares["both"] - F_AMOUNT) < 1e-12 and held["remainder alone"] == 1.0 and held["local alone"] == 0.0
    assert 0.96 < at_tenth < 0.98 and 0.995 < today < 0.998
    out.update(held_exponent=p, held_numeric=numeric, held_readings=held, baryons_at_tenth=at_tenth,
               baryons_today=today)
    report(f"D4 the share 1 - 1/(2 pi) = {F_AMOUNT:.4f}: with the baryons held the amount grows as a^{p:.4f} "
           f"(integrated {numeric:.4f}. The remainder alone gives {held['remainder alone']:.0f} and the local share "
           f"alone {held['local alone']:.0f}). Released "
           f"smooth at a = 1/1100 the baryons reach {at_tenth:.2f} of the amount's excess at a = 1/10 and "
           f"{today:.3f} today")

    report("RESULT: one scale free growth factor, D = a in the matter era and index 0.55 under the constant "
           f"vacuum term. Lambda carries no perturbation. The amount's share sets the held exponent {p:.4f}, and "
           "the baryons join the same mode.")
    return out


if __name__ == "__main__":
    run()
