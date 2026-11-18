"""A Z3 lump in the empty field: clock inside, one cone, weight = energy.

THE READING UNDER TEST. The vacuum is the empty field, the engine's zero
orbit (interference/__init__.py). The clock lives inside the lump, a
self-trapped patch whose phase turns at one rate. Then the vacuum
carries no profile, the carrier's light cone is the only cone, and the
tensor term of composite_metric.py is induced at c_T = c
(carrier_sakharov.py, no condensate). This probe computes the lump.

PREMISES. P1 the carrier of relativistic_carrier.py, three complex
Klein-Gordon fields, hbar = 1. P2 the self-bound sign of
self_bound_vacuum.py, g0 = -1, g1 = 1/sqrt 2, kappa = 1. This sign binds,
with the clock as the phase ground state. P3 the vacuum phi = 0.
P4 for L3's last line, the linearized Einstein response of the composite
metric to the lump's energy momentum tensor.

STATEMENTS (the numbers below check each one).

  L1 (reduction, existence). On the clock orientation phi_a = Phi w^a
     / sqrt 3 the three component force is (2 a R + kappa R^2) psi_a,
     with a = g0/2 - g1/6. So the lump Phi = sigma(r) exp(-i Omega t)
     obeys
        mu f = -(f'' + 2 f'/r)/(2m) + 2 a f^3 + kappa f^5,
     with f = c sqrt(2m) sigma and mu = (Omega^2 - m^2 c^4)/(2 m c^2).
     The profile is the same at every c. Lumps exist for mu0 < mu < 0,
     mu0 = a R0 / 2. For a >= 0 the identity of L3 has one sign and no
     lump exists, so the gas sign (+,+) holds no lump in the empty
     field.
  L2 (clock and charge). With N = int f^2, K = int |grad f|^2/(2m),
     P = int (a f^4 + kappa f^6/3):  Q = Omega N/(m c^2),
     E = (m c^2 + mu) N + K + P. The clock rate Omega lies below m c^2.
     The branch dQ/dOmega < 0 is the stable one. E < m c^2 Q above a
     threshold size. phi = 0 is the global minimum of the potential iff
     m c^2 >= 2 |mu0|.
  L3 (weight). int T_ij d^3x = delta_ij (mu N - P - K/3) = 0 (von
     Laue), at every c. So the active gravitational mass
     int (T_00 + T_kk) is E. Binding energy and clock energy both
     gravitate in full, and gamma - 1 = -4 (mu N - P - K/3)/E.
  L4 (composite metric of the lump). h_00 = 2 Omega^2 sigma^2 and
     h_rr = 2 c^2 sigma'^2. The clock entry is nonzero where the lump is
     and falls as exp(-2 sqrt(-2 m mu) r) outside it.
  L5 (one cone). In 1+1 dimensions, with the full three component
     force, a lump boosted to v moves at v with E = gamma E0 and
     P = gamma E0 v / c^2, for v both below and above the sound speed of
     its own interior. Under band limited noise it stays a clock.

RESULT. Computed on the self bound sign. Z3 lumps with the clock inside
exist in the empty field from N_min = 43.87, and their integrated stress
vanishes to 1e-10 (von Laue), so weight equals energy with binding
included and gamma = 1. In 1+1 dimensions a boosted lump moves at v with
E = gamma E0, also above its interior sound speed, and the gas sign
(+,+) holds no lump in the empty field.

Usage: python3 tests/probes/lump_in_empty_field.py     (about ten seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0, interaction_matrix  # noqa: E402

G0, G1, KAPPA = -1.0, G1_OVER_G0, 1.0
A = G0 / 2 - G1 / 6
R0 = -3 * A / (2 * KAPPA)
MU0 = A * R0 / 2
W3 = np.exp(2j * np.pi * np.arange(3) / 3)


def clock_force_error(rng):
    """L1: M[psi] psi + kappa R^2 psi = (2 a R + kappa R^2) psi on the clock."""
    F = rng.normal() + 1j * rng.normal()
    psi = (F * W3 / math.sqrt(3))[:, None]
    R = abs(F)**2
    force = np.einsum("xab,bx->ax", interaction_matrix(psi, G0, G1), psi) + KAPPA * R**2 * psi
    return np.max(abs(force - (2 * A * R + KAPPA * R**2) * psi))


def _rhs(r, y, mu, a):
    f, fp = y
    return [fp, -2 * fp / r + 2 * (2 * a * f**3 + KAPPA * f**5 - mu * f)]


def _shoot(f0, mu, a, rmax=80.0, dense=False):
    eps = 1e-4
    f2 = 2 * (2 * a * f0**3 + KAPPA * f0**5 - mu * f0) / 3

    def cross(r, y, mu, a):
        return y[0]

    def turn(r, y, mu, a):
        return y[1]
    cross.terminal = turn.terminal = True
    cross.direction, turn.direction = -1, 1
    sol = solve_ivp(_rhs, (eps, rmax), [f0 + 0.5 * f2 * eps**2, f2 * eps], args=(mu, a),
                    events=(cross, turn), rtol=1e-12, atol=1e-14, method="DOP853",
                    dense_output=dense)
    return bool(sol.t_events[0].size), sol


def lump(mu, a=A):
    """Radial lump at chemical potential mu (m = 1): profile and integrals."""
    disc = math.sqrt(a * a + 4 * KAPPA * mu / 3)
    lo = math.sqrt(3 * (-a - disc) / (2 * KAPPA)) * (1 + 1e-9)
    hi = math.sqrt((-a + math.sqrt(a * a + KAPPA * mu)) / KAPPA)
    while True:                                   # bisect on the central value
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        lo, hi = (lo, mid) if _shoot(mid, mu, a)[0] else (mid, hi)
    _, sol = _shoot(lo, mu, a, dense=True)
    r, f, fp = sol.t, sol.y[0], sol.y[1]
    k = math.sqrt(-2 * mu)
    tail = np.where((f > 0) & (f < 1e-4 * lo))[0]
    i = tail[np.argmin(abs(fp[tail] / f[tail] + k + 1 / r[tail]))]
    rc, C = r[i], f[i] * r[i] * math.exp(k * r[i])     # f = C exp(-k r)/r beyond rc
    C2 = (f[i] * r[i])**2

    def profile(x):
        x = np.asarray(x, dtype=float)
        inner = sol.sol(np.clip(x, r[0], rc))[0]
        return np.where(x <= rc, inner, C * np.exp(-k * x) / np.maximum(x, rc))

    def q(g):
        return quad(lambda x: 4 * math.pi * x * x * g(*sol.sol(x)), r[0], rc,
                    limit=400, epsabs=0, epsrel=1e-12)[0]
    N = q(lambda f, fp: f * f) + 4 * math.pi * C2 / (2 * k)
    K = q(lambda f, fp: 0.5 * fp * fp) + 2 * math.pi * C2 * quad(
        lambda x: math.exp(-2 * k * (x - rc)) * (k + 1 / x)**2, rc, np.inf)[0]
    P = q(lambda f, fp: a * f**4 + KAPPA * f**6 / 3)
    wall = float(np.max(fp[:i + 1]**2))
    return dict(mu=mu, f0=lo, N=N, K=K, P=P, E_nr=K + P, s=mu * N - P - K / 3,
                rc=rc, f_rc=f[i], wall=wall, profile=profile)


def lift(d, c=1.0):
    """L2, L3: relativistic charge, energy, weight (m = 1)."""
    mc2 = c * c
    omega = math.sqrt(mc2**2 + 2 * mc2 * d["mu"])
    E = (mc2 + d["mu"]) * d["N"] + d["K"] + d["P"]
    return dict(omega=omega, Q=omega / mc2 * d["N"], E=E,
                binding=E / (omega * d["N"]) - 1, weight=1 + 3 * d["s"] / E,
                gamma_ppn=1 - 4 * d["s"] / E,
                h00=omega**2 * d["f0"]**2 / mc2, hrr=d["wall"],
                h00_far=omega**2 * d["f_rc"]**2 / mc2)


# 1+1 dimensions, m = c = 1: the full three component dynamics.

def profile_1d(mu):
    al, be = A / (2 * mu), math.sqrt(A * A / (4 * mu * mu) + KAPPA / (3 * mu))
    k = 2 * math.sqrt(-2 * mu)

    def f(x):
        return (al + be * np.cosh(k * x))**-0.5

    def fp(x):
        return -0.5 * be * k * np.sinh(k * x) * (al + be * np.cosh(k * x))**-1.5
    return f, fp


def rest_1d(mu):
    f, fp = profile_1d(mu)

    def q(g):
        return quad(g, -60, 60, limit=400, epsabs=0, epsrel=1e-12)[0]
    N, K = q(lambda x: f(x)**2), q(lambda x: 0.5 * fp(x)**2)
    P = q(lambda x: A * f(x)**4 + KAPPA * f(x)**6 / 3)
    return dict(E0=(1 + mu) * N + K + P, virial=mu * N - P + K, omega=math.sqrt(1 + 2 * mu))


def _force(phi):
    psi = math.sqrt(2) * phi
    R = np.sum(abs(psi)**2, axis=0)
    return math.sqrt(2) * (np.einsum("xab,bx->ax", interaction_matrix(psi, G0, G1), psi)
                           + KAPPA * R**2 * psi)


def _measure(phi, pi, x, k):
    dx = x[1] - x[0]
    phix = np.fft.ifft(1j * k * np.fft.fft(phi, axis=1), axis=1)
    psi = math.sqrt(2) * phi
    rho = abs(psi)**2
    R = rho.sum(axis=0)
    w = 0.5 * G0 * R**2 + KAPPA * R**3 / 3 + G1 * sum(
        rho[j] * np.real(np.conj(psi[(j + 1) % 3]) * psi[(j + 2) % 3]) for j in range(3))
    q = -2 * np.sum(np.imag(np.conj(phi) * pi), axis=0)
    modes = np.array([np.sum(np.conj(W3)[:, None]**n * phi, axis=0) / 3 for n in range(3)])
    Q = q.sum() * dx
    return dict(E=(np.sum(abs(pi)**2 + abs(phix)**2 + abs(phi)**2, axis=0) + w).sum() * dx,
                P=-2 * np.sum(np.real(np.conj(pi) * phix)) * dx, Q=Q,
                xc=np.sum(x * q) * dx / Q,
                peak=math.sqrt(np.sum(abs(phi)**2, axis=0).max()),
                off_clock=1 - np.sum(abs(modes[1])**2) / np.sum(abs(modes)**2))


def boosted(mu, v, T=30.0, L=120.0, n=2048, noise=0.0, x0=-15.0):
    """Leapfrog evolution of the lump boosted to v. Samples at nine times."""
    dx = L / n
    x = (np.arange(n) - n // 2) * dx
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    f, fp = profile_1d(mu)
    omega, g = math.sqrt(1 + 2 * mu), 1 / math.sqrt(1 - v * v)
    xi = g * (x - x0)
    phase = np.exp(1j * omega * g * v * (x - x0)) * W3[:, None] / math.sqrt(6)
    phi = f(xi) * phase
    pi = (-g * v * fp(xi) - 1j * omega * g * f(xi)) * phase
    if noise:
        rng = np.random.default_rng(3)
        white = rng.normal(size=phi.shape) + 1j * rng.normal(size=phi.shape)
        smooth = np.fft.ifft(np.fft.fft(white, axis=1) * (abs(k) < 4.0), axis=1)
        phi = phi + noise * np.exp(-(xi / 6)**2) * smooth / abs(smooth).max()
    dt = 0.25 * dx
    steps = 8 * int(round(T / dt / 8))

    def acc(u):
        return np.fft.ifft(-(k * k) * np.fft.fft(u, axis=1), axis=1) - u - _force(u)
    out, a = [dict(t=0.0, **_measure(phi, pi, x, k))], acc(phi)
    for i in range(steps):
        pi = pi + 0.5 * dt * a
        phi = phi + dt * pi
        a = acc(phi)
        pi = pi + 0.5 * dt * a
        if (i + 1) % (steps // 8) == 0:
            out.append(dict(t=(i + 1) * dt, **_measure(phi, pi, x, k)))
    t, xc = np.array([o["t"] for o in out]), np.array([o["xc"] for o in out])
    E, pk = np.array([o["E"] for o in out]), np.array([o["peak"] for o in out])
    return dict(v_fit=np.polyfit(t, xc, 1)[0], gamma=g, E=E[0], P=out[0]["P"],
                dE=abs(E / E[0] - 1).max(), peak=abs(pk / pk[0] - 1).max(),
                off_clock=[o["off_clock"] for o in out])


def run(report=print):
    rng = np.random.default_rng(2)
    force_err = max(clock_force_error(rng) for _ in range(20))
    assert force_err < 1e-12 and A < 0 < R0 and abs(MU0 + 0.286305) < 1e-6
    assert 0.5 - G1 / 6 > 0          # the gas (+,+): mu N - P - K/3 < 0, no lump

    # L1, L2: the family. Minimum of N by a parabola through three points,
    # the E_nr = 0 point by a secant. Then four lumps on the stable branch.
    fr = np.array([0.12, 0.1362, 0.15])
    Ns = np.array([lump(x * MU0)["N"] for x in fr])
    par = np.polyfit(fr, Ns, 2)
    n_min = np.polyval(par, -par[1] / (2 * par[0]))
    ea, eb = lump(0.29 * MU0), lump(0.294 * MU0)
    t = ea["E_nr"] / (ea["E_nr"] - eb["E_nr"])
    n_bound = ea["N"] + t * (eb["N"] - ea["N"])
    assert Ns[1] < min(Ns[0], Ns[2]) and ea["E_nr"] > 0 > eb["E_nr"]
    assert abs(n_min - 43.8675) < 2e-3 and abs(n_bound - 55.673) < 2e-2

    fracs = (0.3, 0.5, 0.7, 0.85)
    lumps = [lump(x * MU0) for x in fracs]
    rel = [lift(d) for d in lumps]
    for d, l in zip(lumps, rel):
        assert abs(d["s"]) < 1e-10 * d["K"]                      # von Laue
        assert abs(l["weight"] - 1) < 1e-10 and abs(l["gamma_ppn"] - 1) < 1e-10
        assert l["omega"] < 1 and l["h00_far"] < 1e-7 * l["h00"]
    Q = [l["Q"] for l in rel]
    assert all(q2 > q1 for q1, q2 in zip(Q, Q[1:]))              # dQ/dOmega < 0
    assert rel[0]["binding"] > 0 > rel[1]["binding"] > rel[2]["binding"] > rel[3]["binding"]
    eotvos = 2 * abs(rel[0]["weight"] - rel[-1]["weight"]) / (rel[0]["weight"] + rel[-1]["weight"])
    passive = [1 / (1 + l["binding"]) for l in rel]              # m Q / E: binding weightless
    eotvos_unbound = 2 * abs(passive[0] - passive[-1]) / (passive[0] + passive[-1])
    assert eotvos < 1e-10 < 0.2 < eotvos_unbound
    c_sound = math.sqrt(-A * R0 / (-A * R0 + 1 + 2 * MU0))       # liquid interior, m = c = 1
    assert 1 > 2 * abs(MU0) and abs(c_sound - 0.75671) < 1e-5

    # L5: boosts in 1+1 dimensions.
    mu1 = 0.6 * MU0
    r1 = rest_1d(mu1)
    assert abs(r1["virial"]) < 1e-10
    runs = {}
    for v in (0.0, 0.6, 0.9):
        b = boosted(mu1, v)
        runs[v] = b
        assert abs(b["v_fit"] - v) < 2e-4 and abs(b["E"] / (b["gamma"] * r1["E0"]) - 1) < 1e-6
        assert b["dE"] < 1e-6 and b["peak"] < 2e-3 and max(b["off_clock"]) < 1e-12
        if v:
            assert abs(b["P"] / (b["gamma"] * r1["E0"] * v) - 1) < 1e-6
    assert 0.9 > c_sound
    noisy = boosted(mu1, 0.0, T=150.0, noise=1e-2, x0=0.0)
    assert max(noisy["off_clock"]) < 2 * noisy["off_clock"][0] and noisy["peak"] < 0.1

    report(f"L1 clock reduction: force error {force_err:.1e}; a = {A:.6f}, R0 = {R0:.6f}, "
           f"mu0 = {MU0:.6f}.")
    report(f"   family: N_min = {n_min:.4f}; E_nr = 0 at N = {n_bound:.3f}.")
    for x, d, l in zip(fracs, lumps, rel):
        report(f"L2 mu/mu0 = {x:.2f}: N = {d['N']:9.3f}, Omega/mc^2 = {l['omega']:.6f}, "
               f"Q = {l['Q']:9.3f}, E/(m c^2 Q) - 1 = {l['binding']:+.5f}; "
               f"L3 (mu N - P - K/3)/K = {d['s'] / d['K']:+.1e}, weight/E - 1 = "
               f"{l['weight'] - 1:+.1e}, gamma - 1 = {l['gamma_ppn'] - 1:+.1e}; "
               f"L4 h_00(0) = {l['h00']:.4f}, max h_rr = {l['hrr']:.4f}, "
               f"h_00(r = {d['rc']:.0f}) = {l['h00_far']:.1e}.")
    report(f"L3 Eotvos parameter, smallest against largest: {eotvos:.1e}. With weightless "
           f"binding it would be {eotvos_unbound:.3f}.")
    report(f"L2 vacuum: phi = 0 is the global minimum for m c^2 >= {2 * abs(MU0):.5f}; "
           f"interior sound speed {c_sound:.5f} c at m = c = 1.")
    for v, b in runs.items():
        mom = f", P/(gamma E0 v) = {b['P'] / (b['gamma'] * r1['E0'] * v):.7f}" if v else ""
        report(f"L5 v = {v:.1f}: v_fit = {b['v_fit']:+.6f}, E/(gamma E0) = "
               f"{b['E'] / (b['gamma'] * r1['E0']):.7f}{mom}, dE/E = {b['dE']:.1e}, "
               f"peak drift {b['peak']:.1e}.")
    report(f"L5 noise 1e-2 over t = 150: off clock share {noisy['off_clock'][0]:.1e} -> "
           f"{noisy['off_clock'][-1]:.1e} (max {max(noisy['off_clock']):.1e}), peak drift "
           f"{noisy['peak']:.1e}.")
    report("RESULT: in the empty field the self-bound sign holds stable Z3 lumps "
           "with the clock inside. Their integrated stress vanishes, so weight equals "
           "energy with binding included and gamma = 1. They move by carrier boosts at "
           "any v < c, including above their own interior sound speed. The gas sign "
           "(+,+) holds no lump in the empty field.")
    return dict(force_error=force_err, n_min=n_min, n_bound=n_bound, lumps=lumps,
                relativistic=rel, eotvos=eotvos, eotvos_unbound=eotvos_unbound,
                c_sound=c_sound, rest_1d=r1, boosts=runs, noisy=noisy)


if __name__ == "__main__":
    run()
