"""The spinning clock lump: a closed circulation of the common phase.

knot_in_lump.py shows the self-bound clock lump holds no pi_3 knot: with
g1 != 0 the only continuous symmetry is the common U(1), and the relative
populations drain. The engine's own knot is a closure, "a single knot
circulating in ONE channel" (interference/words.py). This probe asks
whether the field holds that object: a lump whose common phase winds m
times around a closed loop. Its two Noether charges come from the common
U(1) and from rotations: N and L_z = m N.

PREMISES. The nonrelativistic carrier on the self-bound sign,
g0 = -1, g1 = 1/sqrt 2, kappa = 1, in three dimensions. Its clock
reduction is the cubic quintic equation with a = g0/2 - g1/6 < 0. The
flat top (zero pressure) chemical potential is mu_co = -3 a^2/(4 kappa)
= mu0 (lump_in_empty_field.py).

STATEMENTS (the numbers below check each one).

  S1 (existence). Stationary states psi_a = F(rho, z) e^{i m phi}
     w^a / sqrt 3 found by axisymmetric gradient flow at fixed N: the
     m = 1 torus is bound (E < 0) for N >= 500, with its density peak on
     a ring of radius 4.1, 4.9, 6.1 at N = 500, 1000, 2000 and
     mu/mu_co = 0.59, 0.70, 0.77; m = 2 is bound from N = 500.
  S2 (the stability rule, cited). For the cubic quintic equation in
     three dimensions the m = 1 torus is stable above
     mu/mu_co = 0.69, and no m >= 2 torus is stable. The dangerous
     perturbation splits the ring in two (Mihalache et al., PRL 88,
     073902 (2002); Malomed arXiv:1904.12081 Eq. 19). Multicomponent
     vorticity survives only in the common phase (Kartashov, Malomed,
     Tarruell and Torner, PRA 98, 013612 (2018)).
  S3 (the full three component dynamics). Real time, all three
     components seeded with noise 1e-3 and the ring with 2% n = 2, 3
     dents, absorbing layer, 48^3. RECORD (t = 300, a long run outside the
     suite: evolve(*seeded_3d(torus(N, m)), T=300.0), not asserted):
         m = 1, N = 500:   n = 2 grows to 0.93, the ring splits and
                           sheds 59% of its charge;
         m = 1, N = 1000:  harmonics stay below 0.013, L_z/N = 0.9998,
                           lock 1.0000;
         m = 1, N = 2000:  harmonics below 0.013, L_z/N = 0.9999,
                           lock 1.0000;
         m = 2, N = 2000:  n = 2 grows to 0.55, L_z/N falls to 0.88.
     The rule of S2 holds in the three component field, and the clock
     inside the spinning lump stays locked: the relative sector adds no
     instability.


RESULT. The field holds a closed circulation: one quantum of the common
phase winding around a ring, inside a lump above a threshold charge,
with the clock locked inside it. More winding, or too little charge,
splits it into lumps that carry the charge away. Knotted vortex lines
untie by reconnection into separate rings (Kleckner, Scheeler and
Irvine, Phys. Fluids 26, 091105 (2014)): what persists is the single
closed loop.

Usage: python3 tests/probes/spinning_lump.py     (S1 in under a minute;
run(three_d=True) adds two 3D runs, about ten minutes)
"""

import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import lump_formation as fm  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402

G0, G1, KAPPA = -1.0, lf.G1, 1.0
A = lf.A
MU_CO = -3 * A * A / (4 * KAPPA)
W3 = np.exp(2j * np.pi * np.arange(3) / 3)


def torus(N, m, h=0.2, dtau=0.5, steps=4000, tol=1e-10):
    """S1: axisymmetric stationary state of the clock reduction at fixed N."""
    R0 = (N / (4 / 3 * math.pi * 0.93))**(1 / 3)
    r = h * (np.arange(int((2 * R0 + 12) / h)) + 0.5)
    z = -(R0 + 10) + h * (np.arange(int(2 * (R0 + 10) / h)) + 0.5)
    RR, ZZ = np.meshgrid(r, z, indexing="ij")
    w = (2 * math.pi * RR * h * h).ravel()
    rp, rm = r + h / 2, r - h / 2
    Lr = sp.diags([-(rp + rm) / (r * h * h) - m * m / r**2, rp[:-1] / (r[:-1] * h * h),
                   rm[1:] / (r[1:] * h * h)], [0, 1, -1])
    nz = len(z)
    Lz = sp.diags([-2 * np.ones(nz) / h**2, np.ones(nz - 1) / h**2, np.ones(nz - 1) / h**2],
                  [0, 1, -1])
    L = sp.kron(Lr, sp.eye(nz)) + sp.kron(sp.eye(len(r)), Lz)
    lu = spla.splu((sp.eye(L.shape[0]) - 0.5 * dtau * L).tocsc())
    seed = (RR / (0.6 * R0))**m if m else np.ones_like(RR)
    f = (seed * np.exp(-(np.sqrt(RR**2 + ZZ**2) / (R0 + 1))**4)).ravel()

    def norm(f):
        return f * math.sqrt(N / np.sum(w * f * f))

    def energy(f):
        g = f.reshape(RR.shape)
        fr, fz = np.diff(g, axis=0) / h, np.diff(g, axis=1) / h
        kin = 0.5 * (np.sum(2 * math.pi * (RR[:-1] + h / 2) * h * h * fr**2)
                     + np.sum(2 * math.pi * RR[:, :-1] * h * h * fz**2)
                     + np.sum(w.reshape(RR.shape) * m * m * g**2 / RR**2))
        R = g * g
        return kin + np.sum(w.reshape(RR.shape) * (A * R**2 + KAPPA * R**3 / 3))
    f, last = norm(f), None
    for i in range(steps):
        R = f * f
        f = norm(lu.solve(f - dtau * (2 * A * R + KAPPA * R * R) * f))
        if i % 200 == 0:
            e = energy(f)
            if last is not None and abs(e - last) < tol * abs(e):
                break
            last = e
    R = f * f
    mu = np.sum(w * f * (-0.5 * (L @ f) + (2 * A * R + KAPPA * R * R) * f)) / N
    g = f.reshape(RR.shape)
    i, _ = np.unravel_index(np.argmax(g), g.shape)
    return dict(N=N, m=m, E=energy(f), mu=float(mu), ring=float(r[i]) if m else 0.0,
                r=r, z=z, f=g)


def seeded_3d(t, n=48, L=32.0, noise=1e-3, dent=0.02, seed=0):
    """S3: the torus on a 3D grid, all three components perturbed."""
    dx = L / n
    xs = (np.arange(n) - n // 2) * dx
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    RHO, PHI = np.sqrt(X**2 + Y**2), np.arctan2(Y, X)
    F = RegularGridInterpolator((t["r"], t["z"]), t["f"], bounds_error=False, fill_value=0.0)(
        np.stack([RHO.ravel(), Z.ravel()], -1)).reshape(X.shape)
    F = F * np.exp(1j * t["m"] * PHI) * (1 + dent * np.cos(2 * PHI) + dent * np.cos(3 * PHI))
    rng = np.random.default_rng(seed)
    psi = np.stack([F * w / math.sqrt(3) for w in W3])
    psi = psi + noise * (rng.normal(size=psi.shape) + 1j * rng.normal(size=psi.shape)) * (abs(F) > 0.05)
    psi *= math.sqrt(t["N"] / (np.sum(abs(psi)**2) * dx**3))
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return psi, dict(X=X, Y=Y, Z=Z, RHO=RHO, PHI=PHI, KX=KX, KY=KY, K2=KX**2 + KY**2 + KZ**2,
                     dx=dx, L=L)


def observe(psi, g):
    dx = g["dx"]
    rho = abs(psi)**2
    R = rho.sum(0)
    Lz = 0.0
    for a in range(3):
        P = np.fft.fftn(psi[a])
        dy, dxp = np.fft.ifftn(1j * g["KY"] * P), np.fft.ifftn(1j * g["KX"] * P)
        Lz += float(np.real(np.sum(np.conj(psi[a]) * (-1j) * (g["X"] * dy - g["Y"] * dxp))) * dx**3)
    ring = (abs(g["Z"]) < 2.0) & (g["RHO"] > 1.0)
    wts = R[ring]
    harm = [float(abs(np.sum(wts * np.exp(-1j * k * g["PHI"][ring]))) / np.sum(wts))
            for k in (1, 2, 3, 4)]
    ex = sum(rho[j] * np.real(np.conj(psi[(j + 1) % 3]) * psi[(j + 2) % 3]) for j in range(3))
    return dict(N=float(R.sum() * dx**3), Lz=Lz, harm=harm, lock=float(-6 * ex.sum() / (R * R).sum()))


def evolve(psi, g, T, dt=0.05, every=10.0):
    """lump_formation's integrating factor RK4 in 3D with an absorbing shell."""
    half = np.exp(-0.25j * g["K2"] * dt)
    r = np.sqrt(g["X"]**2 + g["Y"]**2 + g["Z"]**2)
    mask = np.exp(-dt * np.clip((r - 0.4 * g["L"]) / (0.1 * g["L"]), 0, 1)**2)

    def flow(u):
        return np.fft.ifftn(np.fft.fftn(u, axes=(1, 2, 3)) * half, axes=(1, 2, 3))

    def rate(u):
        return -1j * fm.force(u, G0, G1, KAPPA)
    rows, per = [dict(t=0.0, **observe(psi, g))], int(round(every / dt))
    for i in range(1, int(round(T / dt)) + 1):
        k1 = rate(psi)
        a = flow(psi)
        k2 = rate(a + 0.5 * dt * flow(k1))
        k3 = rate(a + 0.5 * dt * k2)
        b = flow(a)
        k4 = rate(b + dt * flow(k3))
        psi = (b + dt / 6 * (flow(flow(k1)) + 2 * flow(k2) + 2 * flow(k3) + k4)) * mask
        if i % per == 0:
            rows.append(dict(t=i * dt, **observe(psi, g)))
    return rows


def run(report=print, three_d=False, T=130.0):
    out = {"tori": {}}
    for m in (1, 2):
        for N in (200, 500, 1000, 2000):
            t = torus(N, m)
            out["tori"][(m, N)] = dict(E=t["E"], mu_ratio=t["mu"] / MU_CO, ring=t["ring"])
            report(f"S1 m = {m}, N = {N:4d}: E = {t['E']:+9.3f}, mu/mu_co = {t['mu'] / MU_CO:.3f}, "
                   f"ring radius {t['ring']:.2f}")
    if three_d:
        for N in (500, 2000):
            t0 = time.time()
            psi, g = seeded_3d(torus(N, 1))
            rows = evolve(psi, g, T)
            last = rows[-1]
            out[("3d", N)] = rows
            report(f"S3 m = 1, N = {N}: t = {last['t']:.0f}: max harmonic "
                   f"{max(last['harm']):.3f}, L_z/N = {last['Lz'] / last['N']:.4f}, "
                   f"N = {last['N']:.1f}, lock {last['lock']:.4f} [{time.time() - t0:.0f} s]")
    report("RESULT: one quantum of common phase winding around a ring persists in a lump above a "
           "threshold charge with the clock locked inside. More winding or less charge splits it.")
    return out


if __name__ == "__main__":
    run(three_d="3d" in sys.argv[1:])
