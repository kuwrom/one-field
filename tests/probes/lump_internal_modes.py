"""The internal modes of the 1D clock lump, and why the lock stalls.

lump_formation.py finds that a lump formed from a random cloud turns
toward the clock but does not reach it within the simulated time (F3).
This probe computes the lump's internal modes and measures what each
one does with energy put into it.

PREMISES. The nonrelativistic carrier on the self-bound sign,
g0 = -1, g1 = 1/sqrt 2, kappa = 1 (lump_formation.py), one space
dimension, the same box, absorbing layer and integrator. The lump is
the closed form clock lump at N = 2.14, the charge left by the
formation runs (mu = -0.2536).

STATEMENTS (the numbers below check each one).

  I1 (not near integrability). In one dimension the three component
     cubic equation with the common term alone is the integrable
     Manakov system. On the lump the g1 exchange and the sextic carry
     0.57 of the interaction energy of the common term (0.66 at
     N = 3). So the dynamics is far from that point, and integrability
     is not what holds the excess.
  I2 (spectrum). The linearization J H / 2 about the lump, written
     spectrally on 400 and 600 points, has no growing mode, only zero
     modes in the clock plane, and three relative internal modes below
     the continuum edge |mu|: omega = 0.0393, 0.1726, 0.2376. For each,
     the first harmonic inside the continuum is the 7th, the 2nd and
     the 2nd.
  I3 (radiation damping). An internal mode at amplitude z loses energy
     as d|z|^2/dt = -Gamma |z|^(2K), where K is its first harmonic in
     the continuum. So |z| ~ t^(-1/(2(K-1))) (Cuccagna and Maeda
     arXiv:2009.00573 Lemma 2.6; in one dimension for K = 2 Buslaev
     and Sulem, Ann. IHP AN 20 (2003) 419, Eq. 1.0.16; for vector
     solitons Pelinovsky and Yang, Stud. Appl. Math. 105 (2000) 245).
     Seeded small, the K = 2 mode loses a third of its lock deficit by
     t = 4000. The K = 7 mode keeps its deficit to 1e-5. So the soft
     relative mode does not radiate.
  I4 (large excess). Seeded at the excess of the formation runs
     (+0.026 over the clock lump), the soft mode is nonlinear and does
     shed energy, down to +0.009 by t = 4000, more slowly in each
     window. The excess drains while it is large and stalls as it
     becomes small.
  I5 (the fingerprint, failed as registered). The registered test was
     a line at twice the soft frequency in the lock signal of formed
     lumps. RECORD (t = 600 to 2400, not rerun here): seeds 1, 2 and 7
     show lines at 0.275, 0.060 and 0.074 against twice the soft
     frequency 0.078, 0.078, 0.053. Their mean late lock is 0.44, 0.86
     and 0.79. Within the simulated time, formed lumps sit outside the
     linear neighbourhood of the clock, and the linear mode picture
     does not describe them.

  I6 (the count of bound modes in 3D). The radial linearization of
     the 3D clock lump in partial waves l = 0 to 3, converged between
     400 and 600 points, binds relative modes below the edge |mu|. Their
     number grows with the lump: distinct frequencies 2, 2, 3, 7
     (states, counting 2l + 1: 2, 2, 5, 21) at mu/mu0 = 0.3, 0.5, 0.7,
     0.85 (N = 57, 121, 460, 3155). At mu0/2 they are 0.0406 and 0.1362
     (l = 0). At 0.7 mu0, 0.0570 and 0.1905 (l = 0) and 0.1309 (l = 1).
     The number of bound modes is set by the lump's size, as for any
     well. It is not three. The l = 0 pair at mu0/2 is rechecked below.
     The full table is a RECORD.

RESULT. Near the clock, energy in the soft relative mode stays there,
because its first harmonic in the continuum is the seventh. A large
excess drains nonlinearly and slows as it shrinks. So conservative
dynamics approaches the clock without completing the lock, in any
dimension where the soft mode sits this deep below the continuum. In 3D
the soft relative mode (l = 0) sits at omega/|mu| = 0.28, and its fourth
harmonic is the first in the continuum. The 3D formation run of
lump_formation.py (RECORD) shows no draining of the excess within its
simulated time.

Usage: python3 tests/probes/lump_internal_modes.py     (about two minutes)
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from arena_rulers import to_real  # noqa: E402
import lump_formation as fm  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402
import lump_stability as ls  # noqa: E402

G0, G1, KAPPA = -1.0, lf.G1, 1.0
W3 = lf.W3
J6 = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
_CLOCK = to_real(W3 / math.sqrt(3))
_ICLOCK = to_real(1j * W3 / math.sqrt(3))
P_CLOCK = np.outer(_CLOCK, _CLOCK) + np.outer(_ICLOCK, _ICLOCK)
N_LUMP = 2.14


def integrability(N):
    """I1: interaction energy off the Manakov point over the common term."""
    mu, _ = ls.lump_1d(N)
    f = lf.profile_1d(mu)[0]

    def q(g):
        return quad(g, -60, 60, limit=400)[0]
    common = q(lambda x: 0.5 * G0 * f(x)**4)
    exchange = q(lambda x: -G1 / 6 * f(x)**4)          # the g1 term on the clock
    sextic = q(lambda x: KAPPA * f(x)**6 / 3)
    return (abs(exchange) + abs(sextic)) / abs(common)


def spectrum(mu, L=40.0, n=600):
    """I2: frequencies, clock plane shares and real space vectors below |mu|."""
    h = 2 * L / n
    x = -L + h * (np.arange(n) + 0.5)
    f = lf.profile_1d(mu)[0](x)
    W = np.array([ls.local_hessian(fx, mu, G1) for fx in f])
    k = 2 * np.pi * np.fft.fftfreq(n, d=h)
    K2 = np.real(np.fft.ifft(np.diag(k**2) @ np.fft.fft(np.eye(n), axis=0), axis=0))
    H = np.kron(K2, np.eye(6))
    for i in range(n):
        H[6 * i:6 * i + 6, 6 * i:6 * i + 6] += W[i]
    ev, V = np.linalg.eig(0.5 * np.kron(np.eye(n), J6) @ H)
    growth = float(np.max(abs(ev.real)))
    modes = []
    for lam, v in zip(ev, V.T):
        if lam.imag < -1e-9 or abs(lam) > abs(mu):
            continue
        u = v.reshape(n, 6)
        share = float(np.sum(np.abs(u @ P_CLOCK)**2) / np.sum(np.abs(u)**2))
        modes.append(dict(omega=float(lam.imag), share=share, x=x, u=u))
    return growth, sorted(modes, key=lambda m: m["omega"])


def bound_3d(frac, ell, n=400, rmax=30.0):
    """I6: converged grid relative bound frequencies of the 3D lump in one partial wave."""
    d = lf.lump(frac * lf.MU0)
    mu = d["mu"]
    h = rmax / (n + 1)
    r = h * np.arange(1, n + 1)
    f = d["profile"](r)
    grid = np.linspace(0.0, f.max() * 1.001, 240)
    table = np.array([ls.local_hessian(x, mu, G1) for x in grid])
    W = np.empty((n, 6, 6))
    for i in range(6):
        for j in range(6):
            W[:, i, j] = np.interp(f, grid, table[:, i, j])
    D = (np.diag(2 * np.ones(n)) - np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / h**2
    H = np.kron(D + np.diag(ell * (ell + 1) / r**2), np.eye(6))
    for i in range(n):
        H[6 * i:6 * i + 6, 6 * i:6 * i + 6] += W[i]
    ev, V = np.linalg.eig(0.5 * np.kron(np.eye(n), J6) @ H)
    out = []
    for lam, v in zip(ev, V.T):
        if 1e-3 < lam.imag < 0.995 * abs(mu):
            u = v.reshape(n, 6)
            if np.sum(np.abs(u @ P_CLOCK)**2) / np.sum(np.abs(u)**2) < 0.05:
                out.append(float(lam.imag))
    merged = []
    for w in sorted(out):
        if not merged or w - merged[-1] > 2e-3:
            merged.append(w)
    return merged


def harmonic(omega, edge):
    return next(k for k in range(1, 64) if k * omega > edge)


def seeded(mu, mode, eps):
    """The clock lump plus eps times the real part of one internal mode."""
    x, _, dx = fm.grid()
    psi = (lf.profile_1d(mu)[0](x) * (W3 / math.sqrt(3))[:, None]).astype(complex)
    re = np.array([np.interp(x, mode["x"], mode["u"][:, j].real) for j in range(6)])
    d = re[:3] + 1j * re[3:]
    return psi + eps * d / math.sqrt(np.sum(abs(d)**2) * dx)


def evolve(psi, T, every=2.0):
    """lump_formation's integrating factor RK4 with its absorbing layer."""
    x, k, dx = fm.grid()
    dt = fm.DT
    half = np.exp(-0.25j * k * k * dt)
    edge = np.clip((abs(x) - 0.36 * fm.BOX) / (0.14 * fm.BOX), 0, 1)
    mask = np.exp(-dt * edge**2)

    def flow(u):
        return np.fft.ifft(np.fft.fft(u, axis=1) * half, axis=1)

    def rate(u):
        return -1j * fm.force(u, G0, G1, KAPPA)
    skip, steps = int(round(every / dt)), int(round(T / dt))
    rows = []
    for i in range(steps + 1):
        if i % skip == 0:
            m = fm.measure(psi, G0, G1, KAPPA)
            rows.append((i * dt, m["lock"], m["E_w"], m["N_w"]))
        if i == steps:
            break
        k1 = rate(psi)
        a = flow(psi)
        k2 = rate(a + 0.5 * dt * flow(k1))
        k3 = rate(a + 0.5 * dt * k2)
        b = flow(a)
        k4 = rate(b + dt * flow(k3))
        psi = (b + dt / 6 * (flow(flow(k1)) + 2 * flow(k2) + 2 * flow(k3) + k4)) * mask
    return np.array(rows)


def windows(rows, edges):
    """Mean lock deficit and mean excess over the clock lump, per window."""
    t, lock, E, N = rows.T
    deficit, excess = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (t >= a) & (t < b)
        deficit.append(float(np.mean(1 - lock[sel])))
        excess.append(float(np.mean(E[sel]) - ls.lump_1d(float(np.mean(N[sel])))[1]))
    return deficit, excess


def run(report=print, T=4000.0):
    out = {"off_manakov": {N: integrability(N) for N in (2.14, 3.0)}}
    mu, _ = ls.lump_1d(N_LUMP)
    edge = abs(mu)
    growth_coarse, coarse = spectrum(mu, n=400)
    growth, modes = spectrum(mu)
    relative = [m for m in modes if m["share"] < 0.05]
    clock_plane = [m for m in modes if m["share"] > 0.95]
    out.update(mu=mu, growth=max(growth, growth_coarse),
               omegas=[m["omega"] for m in relative],
               omegas_coarse=[m["omega"] for m in coarse if m["share"] < 0.05],
               clock_plane_max=max((m["omega"] for m in clock_plane), default=0.0),
               harmonics=[harmonic(m["omega"], edge) for m in relative])
    edges = np.array([0.0, 250.0, 500.0, 1000.0, 2000.0, 3000.0, T])
    runs = {"soft small": (0, 0.12), "hard small": (1, 0.12), "soft large": (0, 0.75)}
    for label, (index, eps) in runs.items():
        rows = evolve(seeded(mu, relative[index], eps), T)
        deficit, excess = windows(rows, edges)
        out[label] = dict(deficit=deficit, excess=excess)
    report(f"I1 interaction energy off the Manakov point / common term: "
           f"{out['off_manakov'][2.14]:.3f} (N = 2.14), {out['off_manakov'][3.0]:.3f} (N = 3).")
    report(f"I2 lump at mu = {mu:.4f}: largest real part {out['growth']:.1e}. Clock plane only "
           f"zero modes (largest {out['clock_plane_max']:.1e}). Relative modes "
           f"{np.round(out['omegas'], 4)} (400 points: {np.round(out['omegas_coarse'], 4)}) "
           f"below the edge {edge:.4f}. First harmonics in the continuum {out['harmonics']}.")
    for label in runs:
        r = out[label]
        report(f"I3/I4 {label:10s}: lock deficit by window {np.round(r['deficit'], 5)}; "
               f"excess over the clock lump {np.round(r['excess'], 5)}.")
    out["bound_3d"] = bound_3d(0.5, 0)
    report(f"I6 3D lump at mu0/2, l = 0: bound relative modes {np.round(out['bound_3d'], 4)} "
           "(RECORD: 2, 2, 3, 7 distinct frequencies at mu/mu0 = 0.3, 0.5, 0.7, 0.85).")
    report("RESULT: the soft relative mode does not radiate (its first harmonic in the continuum "
           "is the seventh). A large excess drains and slows. Within the simulated time, formed "
           "lumps sit outside the linear neighbourhood of the clock (I5, failed fingerprint).")
    return out


if __name__ == "__main__":
    run()
