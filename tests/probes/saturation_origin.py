"""What saturates the self-bound lump, and what the gates take from it.

THE QUESTION. On the self-bound sign the common term of the carrier
binds (a = g0/2 - g1/6 < 0 on the clock) and a saturating term stops
the collapse. The repository carries the total density sextic
(kappa/3) R^3 (self_bound_vacuum.py). Does the carrier supply that
term itself, from its dispersion or from its own fluctuations? If not,
what is its coefficient, and which of its properties do the gates use?

PREMISES. The nonrelativistic limit of the carrier (hbar = m = 1),
the self-bound sign g0 = -1, g1 = 1/sqrt 2 (nls_soliton.G1_OVER_G0),
the uniform clock psi_a = sqrt(R/3) w^a, and for O3 the Gaussian
fluctuations about it (one loop).

STATEMENTS (the numbers below check each one).

  O1 (dispersion). The relativistic lift of the uniform clock
     (relativistic_carrier.py: mu_tot^2 = m^2 c^4 + 2 m c^2 mu_w with
     mu_w = 2 a R) gives mu_NR(R) = m c^2 (sqrt(1 + 4 a R / m c^2) - 1)
     = 2 a R - 2 a^2 R^2 / m c^2 + ...: the term the carrier's
     dispersion adds at the sextic order is attractive. V is U(1)
     invariant, so the envelope expansion generates no harmonic and no
     quintic of its own (relativistic_carrier.py, Sec. 1). The real
     nonlinear Klein-Gordon route to a cubic quintic envelope
     (Sedletsky and Gandzha arXiv:2212.03316) is closed to a complex
     carrier. The cubic only relativistic lump has dQ/dOmega > 0
     (lump_in_empty_field.py).
  O2 (one loop, static). The real Hessian of the quartic on the clock
     has eigenvalues (-7.414, -2.471, -2.307, -2.307, -1.457, -1.457) R,
     all negative, and the one loop sextic of a static background is
     proportional to the sum of their cubes, -453 R^3 / (192 pi^2 m^2)
     per the Coleman-Weinberg expansion: attractive.
  O3 (one loop, zero point). About the uniform clock with no
     saturating term the relative branches are the one_cone.py pair,
     omega^2 = (eps + 3r/4)(eps + 7r/4) up to the shifts -/+ 3r/4,
     r = g1 R/3, and the common branch is omega^2 = eps (eps + 4 a R),
     imaginary for every eps < 4 |a| R: the uniform clock liquid has no
     Gaussian ground state without a classical saturating term. The
     zero point energy of the relative pair is C R^{5/2} with
     C = 8.44e-4 (the Lee-Huang-Yang form, Petrov arXiv:1506.08419
     Eq. 4, validated on the one component value 8/(15 pi^2)). If the
     imaginary common branch is dropped, as Petrov does, the equilibrium
     it sets is at R0 = 2.4e5 with depth R0 xi^3 = 1.5e-3,
     xi = 1/sqrt(2 |a| R0). That is fewer than one quantum per healing
     volume. No classical lump exists there, since a classical lump
     needs R xi^3 >> 1. The fluctuations are not the saturating term. Under
     the node event reading of the proton lifetime the depth is above
     196/c (substrate_depth.py). At that depth their correction to the
     equation of state is below 1e-5 c. If the common mode carries the
     gauge coupling, alpha_G2 = 1/(24 pi) fixes the depth at 3 in the
     range convention xi = 1/(2 c_s) (substrate_coupling.py), which is
     8.49 in this probe's convention xi^2 = 1/(2|a|R0) (channel_
     response.py X3). The correction is then 1.2e-4 and the sextic
     coefficient kappa = 126 (O4).
  O4 (the coefficient is the depth). With the sextic, R0 = 3|a|/(2
     kappa) and the depth is R0 xi^3 = sqrt(kappa) / (2 sqrt 3 a^2):
     kappa enters the lump only through this number. The lump's shape is
     fixed, and psi -> s psi, x -> x/s is the rescaling. Under the
     node event reading the floor is sqrt(kappa)/a^2 > 679/c. The saturating coefficient is the depth
     already registered as premise #1 of the physical layer
     (substrate_depth.py), not a second input.
  O5 (the form). For any saturating function of the total density,
     R^{5/2}, R^3, R^4, the clock lump exists above a zero pressure
     point, obeys mu N - P - K/3 = 0 (the identity behind weight =
     energy, lump_in_empty_field.py), and has the same stable count:
     one negative direction and the phase zero at l = 0, the
     translation zero at l = 1, every relative direction positive. The
     cubic norm sextic 9 kappa |psi_0 psi_1 psi_2|^2, equal to
     (kappa/3) R^3 on the clock, opens the relative sector instead:
     the uniform relative Hessian at R0 has a pair at -2.71 against
     +0.15 for the total density form, at every kappa (both scale as
     1/kappa), and the lump acquires negative pairs at l = 0 and
     l = 1. The clock survives a cubic norm admixture x only below
     x = 0.053. The saturating term acts on the common mode: the
     gates use its existence and its depth, not its power.

RESULT. The carrier's dispersion and its own fluctuations do not
saturate the lump. The saturating term is a common mode term of the
substrate, its coefficient is the depth, and its power is not used.

Usage: python3 tests/probes/saturation_origin.py     (about a minute)
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from arena_rulers import potential, to_complex, to_real  # noqa: E402
from nls_soliton import G1_OVER_G0  # noqa: E402
from self_bound_vacuum import J6, clock, hessian  # noqa: E402

G0, G1 = -1.0, G1_OVER_G0
A = G0 / 2 - G1 / 6
W3 = np.exp(2j * np.pi * np.arange(3) / 3)
FLOOR = 196.0                                   # c R xi^3 > 196, substrate_depth.py

FORMS = {"R^5/2": (lambda R: R**2.5, lambda R: 2.5 * R**1.5),
         "R^3": (lambda R: R**3 / 3, lambda R: R**2),
         "R^4": (lambda R: R**4 / 4, lambda R: R**3)}


def quartic(y):
    return potential(y, G0, G1)


def eos(sat, R):
    """Energy density, chemical potential, pressure of the clock."""
    e, mu = A * R * R + sat[0](R), 2 * A * R + sat[1](R)
    return e, mu, R * mu - e


def zero_pressure(sat):
    R0 = brentq(lambda R: eos(sat, R)[2], 1e-6, 50.0)
    return R0, eos(sat, R0)[1]


def depth(R0):
    return R0 / math.sqrt(2 * abs(A) * R0)**3


def lift_mu(R, mc2):
    """O1: NR chemical potential of the uniform clock from the relativistic shell."""
    return mc2 * (math.sqrt(1 + 4 * A * R / mc2) - 1)


def lhy_integral(l1, l2):
    """int sqrt(s) [sqrt((s+l1)(s+l2)) - s - (l1+l2)/2 + (l1-l2)^2/(8s)] ds."""
    f = lambda s: math.sqrt(s) * (math.sqrt((s + l1) * (s + l2)) - s - (l1 + l2) / 2
                                  + (l1 - l2)**2 / (8 * s))
    return quad(f, 0, np.inf, limit=800)[0]


def branches(R, eps):
    """Bogoliubov frequencies about the uniform clock, no saturating term."""
    mu = 2 * A * R
    H = hessian(lambda y: quartic(y) - mu * float(y @ y), clock(R))
    ev = np.linalg.eigvals(0.5 * J6 @ (2 * eps * np.eye(6) + H)) * 1j
    return np.sort_complex(ev)


def _rhs(r, y, mu, sat):
    f, fp = y
    R = f * f
    return [fp, -2 * fp / r + 2 * (2 * A * R * f + sat[1](R) * f - mu * f)]


def _shoot(f0, mu, sat, rmax=80.0, dense=False):
    eps = 1e-4
    R = f0 * f0
    f2 = 2 * (2 * A * R * f0 + sat[1](R) * f0 - mu * f0) / 3

    def cross(r, y, mu, sat):
        return y[0]

    def turn(r, y, mu, sat):
        return y[1]
    cross.terminal = turn.terminal = True
    cross.direction, turn.direction = -1, 1
    sol = solve_ivp(_rhs, (eps, rmax), [f0 + 0.5 * f2 * eps**2, f2 * eps], args=(mu, sat),
                    events=(cross, turn), rtol=1e-12, atol=1e-14, method="DOP853",
                    dense_output=dense)
    return bool(sol.t_events[0].size), sol


def lump(mu, sat):
    """Radial clock lump for a saturating term sat = (e_sat(R), e_sat'(R))."""
    grid = np.geomspace(1e-6, 50.0, 4000)
    R_min = grid[np.argmin([eos(sat, R)[1] for R in grid])]
    R_hi = brentq(lambda R: eos(sat, R)[1] - mu, R_min, 50.0)
    R_lo = brentq(lambda R: eos(sat, R)[0] - mu * R, 1e-9, R_hi * (1 - 1e-9))
    lo, hi = math.sqrt(R_lo) * (1 + 1e-9), math.sqrt(R_hi)
    while True:
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        lo, hi = (lo, mid) if _shoot(mid, mu, sat)[0] else (mid, hi)
    _, sol = _shoot(lo, mu, sat, dense=True)
    r, f, fp = sol.t, sol.y[0], sol.y[1]
    k = math.sqrt(-2 * mu)
    tail = np.where((f > 0) & (f < 1e-4 * lo))[0]
    i = tail[np.argmin(abs(fp[tail] / f[tail] + k + 1 / r[tail]))]
    rc, C2 = r[i], (f[i] * r[i])**2

    def q(g):
        return quad(lambda x: 4 * math.pi * x * x * g(*sol.sol(x)), r[0], rc,
                    limit=400, epsabs=0, epsrel=1e-12)[0]
    N = q(lambda f, fp: f * f) + 4 * math.pi * C2 / (2 * k)
    K = q(lambda f, fp: 0.5 * fp * fp) + 2 * math.pi * C2 * quad(
        lambda x: math.exp(-2 * k * (x - rc)) * (k + 1 / x)**2, rc, np.inf)[0]
    P = q(lambda f, fp: A * f**4 + sat[0](f * f))

    def profile(x):
        x = np.asarray(x, dtype=float)
        inner = sol.sol(np.clip(x, r[0], rc))[0]
        return np.where(x <= rc, inner, math.sqrt(C2) * np.exp(-k * x) / np.maximum(x, rc))
    return dict(mu=mu, N=N, K=K, P=P, E=K + P, laue=(mu * N - P - K / 3) / K, profile=profile)


def total_density(sat):
    return lambda y: sat[0](float(y @ y))


def cubic_norm(kappa):
    """9 kappa |psi_0 psi_1 psi_2|^2: equal to (kappa/3) R^3 on the clock."""
    return lambda y: 9 * kappa * abs(np.prod(to_complex(y)))**2


def local_hessian(f, mu, sat6):
    x = to_real(f * W3 / math.sqrt(3))
    return hessian(lambda y: quartic(y) + sat6(y) - mu * float(y @ y), x)


def partial_wave(profile, mu, sat6, ell, rmax=30.0, n=900, k=10):
    """Lowest eigenvalues of the second variation in the partial wave l."""
    h = rmax / (n + 1)
    r = h * np.arange(1, n + 1)
    f = profile(r)
    grid = np.linspace(0.0, f.max() * 1.001, 200)
    table = np.array([local_hessian(x, mu, sat6) for x in grid])
    Wr = np.empty((n, 6, 6))
    for i in range(6):
        for j in range(6):
            Wr[:, i, j] = np.interp(f, grid, table[:, i, j])
    N = 6 * n
    main, blocks = np.zeros(N), {}
    for i in range(6):
        main[i::6] = 2 / h**2 + ell * (ell + 1) / r**2 + Wr[:, i, i]
        for j in range(i + 1, 6):
            blocks[j - i] = blocks.get(j - i, np.zeros(N))
            blocks[j - i][i::6] = Wr[:, i, j]
    off = -np.ones(N - 6) / h**2
    data, offsets = [main, off, off], [0, 6, -6]
    for d, v in blocks.items():
        data += [v[:N - d], v[:N - d]]
        offsets += [d, -d]
    L = diags(data, offsets, format="csc")
    return np.sort(eigsh(L, k=k, sigma=-2.0, which="LM", return_eigenvectors=False))


def count(vals, tol=1e-3):
    return int(np.sum(vals < -tol)), int(np.sum(abs(vals) <= tol))


def run(report=print):
    out = {}
    # O1: the dispersion's own sextic order term
    mc2 = 1.0
    h = 1e-3
    kin = (lift_mu(2 * h, mc2) - 2 * lift_mu(h, mc2) + lift_mu(0.0, mc2)) / (2 * h * h)
    out["kinematic_R2"] = kin
    report(f"O1 dispersion: mu_NR = 2 a R + ({kin:+.4f}) R^2 + ..., closed form "
           f"-2 a^2/m c^2 = {-2 * A * A / mc2:+.4f}: attractive.")

    # O2: static one loop
    ev = np.linalg.eigvalsh(hessian(quartic, clock(1.0)))
    out["hessian_eigs"], out["sum_cubes"] = ev, float(np.sum(ev**3))
    report(f"O2 static one loop: quartic Hessian on the clock / R = {np.round(ev, 3)}, "
           f"sum of cubes {out['sum_cubes']:+.1f}: the one loop sextic is attractive.")

    # O3: zero point backreaction
    R = 1.0
    r = G1 * R / 3
    w = branches(R, 1.0)
    closed = sorted([math.sqrt((1 + 5 * r / 4)**2 - r * r / 4) + s * 3 * r / 4 for s in (1, -1)])
    rel = sorted(x.real for x in w if abs(x.imag) < 1e-9 and x.real > 0)
    common = [x for x in w if abs(x.imag) > 1e-9 and x.imag > 0][0]
    out["branch_error"] = float(np.max(abs(np.array(rel) - closed)))
    out["common_imag"] = float(common.imag)
    I_one = lhy_integral(0.0, 2.0)
    out["lhy_check"] = I_one / (32 / (15 * math.sqrt(2)))
    C = math.sqrt(2) / (4 * math.pi**2) * 2 * lhy_integral(3 * r / 4, 7 * r / 4)
    R0_q = (2 * abs(A) / (3 * C))**2
    out["C_lhy"], out["R0_lhy"], out["depth_lhy"] = C, R0_q, depth(R0_q)
    report(f"O3 zero point: relative branches at eps = 1 {np.round(rel, 5)} vs one_cone "
           f"{np.round(closed, 5)}. Common branch {common:.4f} (imaginary). One component "
           f"check {out['lhy_check']:.6f}; C = {C:.3e} R^5/2. The equilibrium it would set: "
           f"R0 = {R0_q:.3e}, depth R0 xi^3 = {out['depth_lhy']:.2e} against the floor "
           f"{FLOOR:.0f}/c.")
    out["correction_at_floor"] = C / (abs(A) * FLOOR * (2 * abs(A))**1.5)
    report(f"   at the floor the zero point correction to the equation of state is "
           f"{out['correction_at_floor']:.1e} c.")

    # O4: the coefficient is the depth
    sat3 = FORMS["R^3"]
    R0, mu0 = zero_pressure(sat3)
    out["depth_formula"] = depth(R0) / (1 / (2 * math.sqrt(3) * A * A))
    out["kappa_floor"] = FLOOR * 2 * math.sqrt(3)
    report(f"O4 depth: R0 = {R0:.4f}, R0 xi^3 = {depth(R0):.4f} = sqrt(kappa)/(2 sqrt3 a^2) "
           f"x {out['depth_formula']:.6f}. The floor reads sqrt(kappa)/a^2 > "
           f"{out['kappa_floor']:.0f}/c.")

    # O5: the form
    out["forms"] = {}
    for name, sat in FORMS.items():
        R0, mu0 = zero_pressure(sat)
        d = lump(0.5 * mu0, sat)
        e0 = partial_wave(d["profile"], d["mu"], total_density(sat), 0)
        e1 = partial_wave(d["profile"], d["mu"], total_density(sat), 1)
        rel_min = float(np.sort(e0)[2])
        out["forms"][name] = dict(R0=R0, laue=d["laue"], l0=count(e0), l1=count(e1),
                                  relative_min=rel_min)
        report(f"O5 {name}: R0 = {R0:.4f}, lump N = {d['N']:.2f}, (mu N - P - K/3)/K = "
               f"{d['laue']:.1e}, counts l=0 {count(e0)}, l=1 {count(e1)}, lowest relative "
               f"{rel_min:+.4f}")
    R0, mu0 = zero_pressure(sat3)
    d = lump(0.5 * mu0, sat3)
    e0 = partial_wave(d["profile"], d["mu"], cubic_norm(1.0), 0)
    e1 = partial_wave(d["profile"], d["mu"], cubic_norm(1.0), 1)
    out["norm_l0"], out["norm_l1"] = count(e0), count(e1)
    tot = np.linalg.eigvalsh(local_hessian(math.sqrt(R0), mu0, total_density(sat3)))
    nrm = np.linalg.eigvalsh(local_hessian(math.sqrt(R0), mu0, cubic_norm(1.0)))
    out["uniform_total"], out["uniform_norm"] = float(tot[1]), float(nrm[0])
    # admixture x of the cubic norm: relative eigenvalue linear in x
    out["x_critical"] = tot[1] / (tot[1] - nrm[0])
    for kappa in (0.25, 4.0):
        R0k, mu0k = zero_pressure((lambda R, k=kappa: k * R**3 / 3, lambda R, k=kappa: k * R**2))
        nk = np.linalg.eigvalsh(local_hessian(math.sqrt(R0k), mu0k, cubic_norm(kappa)))
        out[f"uniform_norm_kappa_{kappa}"] = float(nk[0] * kappa)
    report(f"O5 cubic norm: uniform relative Hessian at R0 {nrm[0]:+.4f} (pair) against "
           f"{tot[1]:+.4f}. Kappa x lowest = {out['uniform_norm_kappa_0.25']:+.4f}, "
           f"{out['uniform_norm_kappa_4.0']:+.4f} at kappa = 1/4, 4. Lump counts l=0 "
           f"{count(e0)}, l=1 {count(e1)}. The clock survives an admixture below x = "
           f"{out['x_critical']:.3f}.")
    report("RESULT: the carrier's dispersion and its own fluctuations do not saturate the "
           "lump. The saturating term is a common mode term of the substrate whose "
           "coefficient is the depth and whose power the gates do not use.")
    return out


if __name__ == "__main__":
    run()
