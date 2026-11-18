"""Self-bound vacuum: does the substrate sit at a zero pressure point?

THE GATE. Nothing outside holds the field, so a state that persists
must hold itself together: p(n0) = 0 at a finite density with no external pressure.  At
that point the thermodynamic potential vanishes, eps - mu n = -p = 0
(Volovik, gr-qc/0604062), the relation dark_sector.py and gravity.py
use inside a lump.  This probe computes which sign of the couplings is
self-bound.

THE GAS SIGN (+,+) (the functional of tests/probes/arena_rulers.py,
alpha = 0, g1 = g0/sqrt 2 > 0):
V = (g0/2) R^2 + g1 sum_k rho_k Re(psi*_{k+1} psi_{k+2}),
R = rho_tot.  On the clock phi_a = sqrt(R/3) w^a the exchange term is
-R^2/6, so eps(R) = a R^2 with a = g0/2 - g1/6, mu = 2 a R and
p = R mu - eps = a R^2 = eps.  With a > 0 the functional is a gas at
every density: it exists only under external pressure, and its
thermodynamic potential is -p.

THE SELF-BOUND SIGN (-,+).  A stable zero pressure point of
eps(R) = a R^2 + b R^3 needs a < 0 < b (p = a R^2 + 2 b R^3 vanishes at
R0 = -a/(2b), and dp/dR = -a R0 > 0 there).  The phase lock is the
g1 term alone: the clock is the ground state at fixed R iff g1 > 0.
Keeping |g1/g0| = 1/sqrt 2, the three sign choices give
    (+,+)  a > 0: a gas;
    (-,-)  a < 0 but the phase lock flips and the clock is not the
           ground state;
    (-,+)  a < 0 with the clock as ground state: with the total density
           sextic (kappa/3) R^3 the substrate is a self-bound liquid.
           This is the committed sign.
The sextic and the quartic in R touch only the common mode: the
relative (matter) Hessian depends on g1 rho alone, so the relative
gaps sqrt(21) r/4 -/+ 3r/4 (one_cone.py) are unchanged.  At R0 the
clock rate mu0 = eps/R0 is fixed by the functional: the density and
the clock are set by the zero pressure point, not by a container.

WHAT IS MEASURED: the gas sign's equation of state. The sign table. For
the liquid (g0 = -1, g1 = 1/sqrt 2, kappa = 1) the zero pressure point,
the clock as phase ground state there, a real Bogoliubov spectrum on a
momentum grid, the relative gaps against the closed form, the phonon
speed against sqrt(-a R0), and the cone ratio c_rel/c_ph.

INTERIOR SPEEDS. c_rel^2 = 5 g1 rho/4.  Gas: c_ph^2 = mu, equal to
c_rel^2 at g1/g0 = 4/3.  Liquid: c_ph^2 = -a R0, equal to c_rel^2 at
g1 = 2|g0|.  At |g1/g0| = 1/sqrt 2 the ratio c_rel/c_ph is 0.621 (gas)
and 0.691 (liquid).  These are internal speeds of the medium.  Gate
C's cone is set by the carrier (tests/probes/lump_in_empty_field.py).

RESULT. On the gas sign (+,+) the clock has p = eps = a R^2 with a > 0,
so it holds nothing without external pressure. At |g1/g0| = 1/sqrt 2
only g0 < 0 < g1 keeps the clock as ground state with a < 0, and with
the total density sextic that sign reaches a zero pressure point where
eps - mu n = 0 to 1e-8. There the Bogoliubov spectrum is real, the
relative gaps keep their closed form, and c_rel/c_ph = 0.691 (0.621 on
the gas sign).

Usage: python3 tests/probes/self_bound_vacuum.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from arena_rulers import potential, to_real  # noqa: E402

W = np.exp(2j * np.pi / 3)
S2 = 1.0 / math.sqrt(2.0)
J6 = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])


def energy(x, g0, g1, kappa):
    return potential(x, g0=g0, g1=g1) + kappa / 3.0 * float(x @ x)**3


def clock(R):
    return to_real(math.sqrt(R / 3.0) * W**np.arange(3))


def eos(g0, g1, kappa, R, h=1e-6):
    """Energy density, chemical potential and pressure of the clock."""
    f = lambda r: energy(clock(r), g0, g1, kappa)
    eps = f(R)
    mu = (f(R + h) - f(R - h)) / (2 * h)
    return eps, mu, R * mu - eps


def clock_is_ground(g0, g1, kappa, R, trials=3000, seed=0):
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(trials):
        w = rng.dirichlet(np.ones(3))
        psi = np.sqrt(w * R) * np.exp(1j * rng.uniform(0, 2 * np.pi, 3))
        best = min(best, energy(to_real(psi), g0, g1, kappa))
    return energy(clock(R), g0, g1, kappa) <= best + 1e-9


def hessian(f, x, h=1e-4):
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei, ej = np.eye(n)[i] * h, np.eye(n)[j] * h
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej)
                       - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return H


def bogoliubov(g0, g1, kappa, R, ps):
    """Frequencies of i dt psi = dE/dpsi* about the clock, per momentum;
    None where an eigenvalue has a real part (instability)."""
    _, mu, _ = eos(g0, g1, kappa, R)
    H = hessian(lambda y: energy(y, g0, g1, kappa) - mu * float(y @ y),
                clock(R))
    out = []
    for p in ps:
        ev = np.linalg.eigvals(0.5 * J6 @ (p * p * np.eye(6) + H))
        out.append(None if np.max(abs(ev.real)) > 1e-6
                   else np.sort(abs(ev.imag))[::2])
    return out


def phonon_c2(g0, g1, kappa, R):
    ps = np.array([0.05, 0.08, 0.12])
    w0 = np.array([min(s) for s in bogoliubov(g0, g1, kappa, R, ps)])
    return float(np.mean((w0**2 - ps**4 / 4) / ps**2))


def run(report=print):
    # gas sign equation of state: a gas at every density
    a_c = 0.5 - S2 / 6
    gas = []
    for R in (0.5, 1.0, 3.0, 10.0):
        eps, mu, p = eos(1.0, S2, 0.0, R)
        assert abs(eps - a_c * R**2) < 1e-9 and abs(p - eps) < 1e-6
        gas.append((R, eps, mu, p))
    assert a_c > 0

    # sign table at |g1/g0| = 1/sqrt 2
    signs = {}
    for g0, g1 in ((1.0, S2), (-1.0, -S2), (-1.0, S2)):
        signs[(g0, g1)] = (g0 / 2 - g1 / 6, clock_is_ground(g0, g1, 0.0, 1.0))
    assert signs[(1.0, S2)][0] > 0 and signs[(1.0, S2)][1]
    assert signs[(-1.0, -S2)][0] < 0 and not signs[(-1.0, -S2)][1]
    assert signs[(-1.0, S2)][0] < 0 and signs[(-1.0, S2)][1]

    # the self-bound liquid
    g0, g1, kappa = -1.0, S2, 1.0
    a = g0 / 2 - g1 / 6
    R0 = -3.0 * a / (2.0 * kappa)
    eps0, mu0, p0 = eos(g0, g1, kappa, R0)
    assert abs(p0) < 1e-8 and abs(eps0 - mu0 * R0) < 1e-8
    ground = clock_is_ground(g0, g1, kappa, R0)
    assert ground
    grid = np.array([0.0, 0.2, 0.5, 1.0, 2.0, 4.0])
    spec = bogoliubov(g0, g1, kappa, R0, grid)
    stable = all(s is not None for s in spec)
    assert stable
    r = g1 * R0 / 3
    gaps = np.sort(spec[0])[1:]
    gaps_expected = np.array([math.sqrt(21) * r / 4 - 0.75 * r,
                              math.sqrt(21) * r / 4 + 0.75 * r])
    assert np.max(abs(gaps - gaps_expected)) < 1e-4
    c2_liquid = phonon_c2(g0, g1, kappa, R0)
    assert abs(c2_liquid - (-a * R0)) < 1e-3
    c_rel2 = 5 * r / 4
    ratio_liquid = math.sqrt(c_rel2 / c2_liquid)
    ratio_gas = math.sqrt((5 * S2 / 4) / (3 - S2))
    assert abs(ratio_liquid - math.sqrt(5 * S2 / (2 * S2 + 6))) < 1e-3
    assert abs(ratio_gas - 0.621) < 1e-3

    for R, eps, mu, p in gas:
        report(f"gas sign, clock R={R:5.1f}: eps={eps:9.4f} mu={mu:8.4f} "
               f"p={p:+9.4f}  (eps - mu n = {-p:+.4f})")
    report(f"gas sign (+,+): a = g0/2 - g1/6 = {a_c:.4f} > 0, p = eps at every "
           "density: a gas, not self-bound.")
    for (s0, s1), (aa, gr) in signs.items():
        report(f"  signs g0={s0:+.0f}, g1={s1:+.3f}: a={aa:+.4f}, clock is "
               f"the phase ground state: {gr}")
    report(f"liquid (g0=-1, g1=+1/sqrt2, kappa=1): R0={R0:.5f}, p={p0:+.1e}, "
           f"eps - mu n = {eps0 - mu0 * R0:+.1e}, clock rate mu0 = {mu0:.5f}")
    report(f"  Bogoliubov spectrum real on the grid: {stable}. Relative gaps "
           f"{np.round(gaps, 5)} = closed form {np.round(gaps_expected, 5)}")
    report(f"  phonon c^2 = {c2_liquid:.5f} = -a R0 = {-a * R0:.5f}; "
           f"c_rel/c_ph = {ratio_liquid:.4f} (gas {ratio_gas:.4f})")
    report("RESULT: the gas sign (+,+) has p = eps > 0 and holds nothing "
           "without external pressure. The sign g0 < 0 < g1 with the "
           "total density sextic is a self-bound liquid that keeps the "
           "clock and the matter gaps and sets eps - mu n = 0 exactly. "
           "c_rel/c_ph is an internal ratio of the medium.")
    return dict(gas=gas, a_gas=a_c, signs=signs, R0=R0,
                p0=p0, interior_term=eps0 - mu0 * R0, mu0=mu0, ground=ground,
                stable=stable, gaps=gaps, gaps_expected=gaps_expected,
                c2_liquid=c2_liquid, ratio_liquid=ratio_liquid,
                ratio_gas=ratio_gas)


if __name__ == "__main__":
    run()
